Project File
~~~~~~~~~~~~

.. list-table:: Table: .prj File lines and their meanings
    :widths: 20 80
    :header-rows: 1

    * - Line
      - Description
    * - ``I,fill2.png``
      - The image file associated with this *.prj* file.
    * - ``D,dim,2``
      - | **Choose either 2D or 2.5D.**
        | 2.5D is the 2D image extruded in the y-direction.
    * - ``D,nx,240``
      - Read from the image file. Do not change.
    * - ``D,ny,1`` 
      - **Number of pixels in the extruded direction if using 2.5D.**
    * - ``D,nz,50``
      - Read from the image file. Do not change.
    * - ``D,dx,1``
      - **Number of meters each pixel represents in the x direction.**
    * - ``D,dy,1``
      - **Number of meters each pixel represents in the y direction.**
    * - ``D,dz,1``
      - **Number of meters each pixel represents in the z direction.**
    * - ``D,cpml,20``
      - **Thickness of absorbing boundary layer. A typical value is 20.**
    * - ``D,nmats,3``
      - Read from the image file. Do not change.
    * - ``D,tfile,``
      - An attenuation processing value that is not yet implemented.
    * - | ``M,1,ice1h,98/197/178,``
        | ``-10,2,910,0,0,TRUE,``
        | ``test.ang``
      - | One comma-separated-values line per material
        | (per color in the model image).
        | **User should change/add the material name (see list),**
        | **temperature (in degrees Celsius), density, porosity, water content,**
        | **whether the material is anisotropic (TRUE or FALSE),**
        | **and if anistropic, the name of the anisotropy file.**
        | **Use a dummy value of 2 for attenuation, recognizing that**
        | **attenuation is not yet incorporated in the calculations.**
        | User should not change the material ID or R/G/B values.
        | Note: Since large density gradients cause numerical instabilities,
        | the density for air must be increased.
        | A value of 400.0 works until a better formulation
        | of the air-rock interface is implemented.
    * - ``S,dt,``
      - Timestep will be calculated automatically.
    * - ``S,time_steps,500``
      - **Decide how many timesteps you want the model to run.**
    * - ``S,x,100``
      - **x-coordinate of the seismic source**
    * - ``S,y,0``
      - **y-coordinate of the seismic source**
    * - ``S,z,0``
      - **z-coordinate of the seismic source**
    * - ``S,f0,60``
      - **Frequency of the seismic source**
    * - ``S,theta,0``
      - **Inclination of the seismic source (+ is down)**
    * - ``S,phi,0``
      - | **Angle of seismic source from x-axis in the x-y plane**
        | **(+ is counterclockwise when viewed from above)**
    * - ``C,0.0,``
      - | Stiffness tensor for each material. User *can*
        | enter or change this manually if desired. If blank,
        | calculated from materials information in the earlier section.
    * - ``E,dt,``
      - Timestep will be calculated automatically.
    * - ``E,time_steps,500``
      - **Number of timesteps to run the model.**
    * - ``E,x,100``
      - **x-coordinate of the radar source**
    * - ``E,y,0``
      - **y-coordinate of the radar source**
    * - ``E,z,0``
      - **z-coordinate of the radar source**
    * - ``E,f0,1e8``
      - | **Frequency of the radar source.**
        | **10-100MHz is a good range to start.**
    * - ``E,theta,0``
      - **Inclination of the radar source (+ is down)**
    * - ``E,phi,0``
      - | **Angle of radar source from x-axis in the x-y plane**
        | **(+ is counterclockwise when viewed from above)**
    * - ``P,0.0,``
      - | Permittivity tensor for each material.
        | User *can* enter or change this manually if desired.
        | If blank, calculated from materials information in the
        | earlier section.
      

* *orientation file*

    A delimited file of one entry of Bunge notation Euler angles per line.
    A typical number of entries is 500 to ensure a smooth data field.