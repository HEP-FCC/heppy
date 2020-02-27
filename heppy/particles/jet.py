import math
import copy
from p4 import P4

def group_pdgid(ptc):
    '''If ptc is a hadron, the pdgid is set to either 211 if it's charged,
    and 130 if it's neutral. Otherwise, the pdgid is the one returned by ptc.pdgid().
    
    @return: the pdgid
    '''
    pdgid = abs(ptc.pdgid())
    if pdgid>100:
        if ptc.q():
            return 211
        else:
            return 130
    else:
        return pdgid

class JetComponent(object):
    '''L{Jet} constituent particle information.
    
    In addition to the values returned by the various methods of this class,
    this class is a list storing each particle added using the L{append} method
    '''
    
    def __init__(self, pdgid):
        '''Create a component for particles of type pdgid.'''
        super(JetComponent, self).__init__()
        self._e = 0
        self._pt = 0
        self._num = 0
        self._pdgid = pdgid
        self._q = 0
        self._particles = list()

    def pdgid(self):
        '''@return: the pdgid'''
        return self._pdgid
        
    def e(self):
        '''@return: the total energy of all particles with this pdgid'''
        return self._e

    def pt(self):
        '''@return: the total pt of all particles with this pdgid'''
        return self._pt

    def num(self):
        '''@return: the number of particles with this pdgid'''
        return self._num
    
    def q(self):
        '''@return: total charge of the particles with this pdgid'''
        return self._q
    
    def particles(self):
        '''@return list of particles'''
        return self._particles
    
    def append(self, ptc):
        '''Append a new particle, incrementing all quantities'''
        pdgid = group_pdgid(ptc)
        if self._pdgid is None:
            self._pdgid = pdgid
        elif pdgid!=self._pdgid:
            raise ValueError('cannot add particles of different type to a component')
        self._particles.append(ptc)
        self._e += ptc.e()
        self._pt += ptc.pt()
        self._q += ptc.q()
        self._num += 1

    def sort(self, *args, **kwargs):
        self._particles.sort(*args, **kwargs)

    def __str__(self):
        header = '\t\tpdgid={pdgid}, n={num:d}, e={e:3.1f}, pt={pt:3.1f}'.format(
            pdgid = self.pdgid(),
            num = self.num(),
            e = self.e(),
            pt = self.pt()
        )
        result = [header]
        if hasattr(self, '_particles'):
            # not there if the component was deepcopied, see below
            ptcs = []
            for ptc in self._particles:
                ptcs.append('\t\t\t{particle}'.format(particle=str(ptc)))
            result.extend(ptcs)
        return '\n'.join(result)
    
    def __deepcopy__(self, memodict={}):
        newone = type(self).__new__(type(self))
        for attr, val in self.__dict__.iteritems():
            if attr not in ['_particles']:
                setattr(newone, attr, copy.deepcopy(val, memodict))
        return newone
        
 
class JetConstituents(object):
    '''Dictionary of constituents.
    
    The dictionary is indexed by the following integer keys:
     - 211: charged hadrons
     - 22: photons
     - 130: neutral hadrons
     - 11: eletrons
     - 13: muons
     
    Each key corresponds to a L{JetComponent}
    '''

    def __init__(self):
        super(JetConstituents, self).__init__()
        all_pdgids = [211, 22, 130, 11, 13, 
                      1, 2 #HF had and em 
                      ]
        self._components = dict()
        for pdgid in all_pdgids:
            self._components[pdgid] = JetComponent(pdgid)
        self.particles = []

    def validate(self, jet_energy, tolerance = 1e-2):
        '''Calls pdb if total component energy != jet energy'''
        tote = sum([comp.e() for comp in self._components.values()]) 
        if abs(jet_energy-tote)>tolerance: 
            import pdb; pdb.set_trace()
    
    def append(self, ptc):
        '''Appends a particle to the list of constituents.'''
        pdgid = group_pdgid(ptc)
        try:
            self._components[pdgid].append(ptc)
        except KeyError:
            msg = '''Particle
            {ptc}
            cannot be added to the jet as it cannot be interpreted
            as a charged hadron, a neutral hadron, a photon, an electron or a muon.
            are you sure that you are using reconstructed particles or
            stable and visible particles from a generator? 
            '''.format(ptc=str(ptc))
            raise ValueError(msg)
        self.particles.append(ptc)
            
    def n_particles(self):
        return sum(comp.num() for comp in self._components.values())
    
    def n_charged_hadrons(self):
        return self._components[211].num()
            
    def __getattr__(self, attr):
        return getattr(self._components, attr)
    
    def __getitem__(self, item):
        return self._components[item]
            
    def sort(self):
        '''Sort constituent particles by decreasing energy.'''
        for ptcs in self._components.values():
            ptcs.sort(key = lambda ptc: ptc.e(), reverse=True)

    def __str__(self):
        return '\n'.join(map(str, self._components.values()))

    def __deepcopy__(self, memodict={}):
        newone = type(self).__new__(type(self))
        for attr, val in self.__dict__.iteritems():
            if attr not in ['particles']:
                setattr(newone, attr, copy.deepcopy(val, memodict))
        return newone

class JetTags(dict):
    '''Dictionary of tags attached to a jet.
    
    The key is the name of the tag, and the value the one of the tag, e.g.:
     - key = 'btag'
     - value = '0.65'  (in case a b tagging algorithm returning a float was used)
     
    Any kind of type can be used for key and value.
    '''

    def summary(self):
        '''@return summary of all tags.'''
        tagstrs = []
        for name, val in sorted(self.iteritems()):
            valstr = '..'
            if hasattr(val, 'summary'):
                valstr = val.summary()
            elif isinstance(val, int):
                valstr = '{val:d}'.format(val=val)
            else:
                try:
                    valstr = '{val:2.1f}'.format(val=val)
                except:
                    pass
            tagstr = '{name}:{val}'.format(name=name, val=valstr)
            tagstrs.append(tagstr)
        return ', '.join(tagstrs)
            
    
class Jet(P4):
    '''Generic jet class.
    
    Attributes:
     - constituents: see L{JetConstituents}
     - tags: see L{JetTags}
    '''

    def __init__(self, *args, **kwargs):
        super(Jet, self).__init__(*args, **kwargs)
        self.constituents = {}
        self.tags = JetTags()
        self._q = None

    def pdgid(self):
        '''needed to behave as a particle.
    
        @return 0
        '''
        return 0

    def q(self):
        if self._q is None:
            self._q = sum(const.q() for const in self.constituents.values())
        return self._q
    
    def __str__(self):
        tmp = '{className} : {p4}, tags={tags}'
        return tmp.format(
            className = self.__class__.__name__,
            p4 = super(Jet, self).__str__(),
            tags = self.tags.summary()
            )
    
    def __repr__(self):
        return str(self)
