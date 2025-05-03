Getting Started
==============================

.. toctree::
   :titlesonly:
   :numbered:
   :hidden:

   installation
   quickstart
   WHATSNEW

``Pryngles`` is a research package, designed and released at the pace of research. For the same reason, it has evolved in a different way compared to professional software.

Due to historical (the history of development) reasons, ``Pryngles`` has two interfaces. These interfaces differ in many aspects; some of them are superficial, others are more profound and are related to the design of the package, even with the physics involved. The interfaces are:

- **System Interface**: This is the newest interface. In this interface, we use the ``snake_case`` naming convention for methods. The interface is loaded using:

  .. code-block:: python

     import pryngles as pr
     sys = pr.System()

- **RingedPlanet Interface**: This was the first interface developed in Pryngles. For the same reason, it is the most tested interface and the one that has been currently used to produce research results. This interface is loaded when the package is imported as:

   .. code-block:: python

     import pryngles as pr
     RP = pr.RingedPlanet()

These interfaces are, in the present version of the package, somewhat compatible. The following commands show how to define a RingedPlanet starting with the System interface:

   .. code-block:: python

      sys = pr.System()

      S = sys.add(kind="Star",
                  physics=dict(radius=Consts.rsun/sys.ul),
                  optics=dict(limb_coeffs=[0.65]))

      P = sys.add(kind="Planet", primary=S,
                  orbit=dict(a=0.2, e=0.0),
                  physics=dict(radius=Consts.rsaturn/sys.ul))

      R = sys.add(kind="Ring", primary=P,
                  physics=dict(fi=1.5, fe=2.5, i=30*Consts.deg))

      O = sys.add(kind="Observer",
                  optics=dict(lamb=90*Consts.deg, beta=90*Consts.deg))

      RP = sys.ensamble_system()
