import pryngles as pr
from pryngles import Consts
sys=pr.System(units=['m', 'kg', 's'])

nspangles=100
S=sys.add(kind="Star",
                    physics=dict(radius=Consts.rsun),
                    optics=dict(limb_coeffs=[0]),
                    orbit=dict(m=1*Consts.msun)
          ,nspangles=100
                   )

O=sys.add(kind="Observer",
                    optics=dict(lamb=0*Consts.deg,beta=0*Consts.deg)
          ,nspangles=100
                   )
P=sys.add(kind="Planet",primary=S,
                    orbit=dict(a=Consts.au,e=0.0, m=1*Consts.mearth),
                    physics=dict(radius=Consts.rearth)
          ,nspangles=100
                   )
R=sys.add(kind="Ring",primary=P,
                    physics=dict(fi=1,fe=1,i=90*Consts.deg)
          ,nspangles=100
                   )

RP = sys.ensamble_system()

print(RP.T)

sys=pr.System(units=['au', 'kg', 's'])

S=sys.add(kind="Star",
                    physics=dict(radius=Consts.rsun/sys.ul),
                    optics=dict(limb_coeffs=[0]),
                    orbit=dict(m=1*Consts.msun)
          ,nspangles=100
                   )

O=sys.add(kind="Observer",
                    optics=dict(lamb=0*Consts.deg,beta=0*Consts.deg)
                   )
P=sys.add(kind="Planet",primary=S,
                    orbit=dict(a=1,e=0.0, m=Consts.mearth),
                    physics=dict(radius=Consts.rearth/sys.ul)
          ,nspangles=100
                   )
R=sys.add(kind="Ring",primary=P,
                    physics=dict(fi=1,fe=1,i=90*Consts.deg)
          ,nspangles=100
                   )

RP = sys.ensamble_system()
print(RP.T)

sys=pr.System(units=['au', 'msun', 'yr'])

S=sys.add(kind="Star",
                    physics=dict(radius=Consts.rsun/sys.ul),
                    optics=dict(limb_coeffs=[0]),
                    orbit=dict(m=Consts.msun/sys.um)
          ,nspangles=100
                   )

O=sys.add(kind="Observer",
                    optics=dict(lamb=0*Consts.deg,beta=0*Consts.deg)
                   )
P=sys.add(kind="Planet",primary=S,
                    orbit=dict(a=Consts.au/sys.ul,e=0.0, m=Consts.mearth),
                    physics=dict(radius=Consts.rearth/sys.ul)
          ,nspangles=100
                   )
R=sys.add(kind="Ring",primary=P,
                    physics=dict(fi=1,fe=1,i=90*Consts.deg)
          ,nspangles=100
                   )

RP = sys.ensamble_system()
print(RP.T)
