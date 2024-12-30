#------------------------------------------------------------------------------------------------#

   _____ ______ _  _   ______                    
  / ____|  ____| || | |  ____|                   
 | (___ | |__  | || |_| |__ ___   __ _ _ __ ___  
  \___ \|  __| |__   _|  __/ _ \ / _` | '_ ` _ \ 
  ____) | |____   | | | | | (_) | (_| | | | | | |
 |_____/|______|  |_| |_|  \___/ \__,_|_| |_| |_|
                                                 
#------------------------------------------------------------------------------------------------#

## Authors 
[                 Ehsan Golab, Sharif University of Technology.                                  ]
[                           Ehsan1996Golab@gmail.com                                             ]                

[                 Nima Shirani, Islamic Azad University.                                         ]
[                           nima.shirani@gmail.com                                               ] 

## Description

The SE4Foam toolbox is a comprehensive simulation suite built upon the robust [**FoamExtend 4.1**] framework. It is designed to address the complex and diverse needs of solar energy research and engineering. SE4Foam extends the capabilities of conventional simulation tools, allowing for the detailed modeling and analysis of a broad spectrum of solar energy systems, including but not limited to:

•	Parabolic trough collectors (PTC)
•	Porous- Parabolic trough collectors (PorousPTC)
•	Photovoltaic (PV) panels
•	Direct absorption collectors
•	Solar water heaters
•	Solar air heaters
•	Solar Chimney

## Key Features
•	[**Versatile Simulation**]: From basic to complex systems, SE4Foam handles a wide range of solar energy applications.

•	[**Advanced Modeling**]: Incorporate modifications in boundary conditions and geometry for a tailored simulation experience.

•	[**Open Source**]: Built on FoamExtend 4.1, SE4Foam promotes collaborative development and innovation.

•	[**Integration of CoolProp**]: Seamlessly integrates the CoolProp library to deliver precise material thermal property computations, a cornerstone for solar energy simulations.

•	[**Adaptive Temperature-Responsive Properties**]: Recognizes the significant temperature variations in solar apparatus, with many physical properties showing a dependency on temperature that surpasses OpenFOAM's standard offerings.

•	[**Innovative Wall Function for Thermal Turbulence**]: Introduces a specialized wall function for turbulent thermal diffusivity, heightening the fidelity of simulations by accurately depicting the heat transfer dynamics at the boundary layer.

•	[**Custom-Designed Boundary Conditions**]: Comes with custom-configured boundary conditions for velocity, heat flux, and temperature, which are essential for resolving issues associated with solar devices.

•	[**Superior Post-Processing Capabilities**]: Outfitted with a broad array of parameters for in-depth post-processing, offering a more comprehensive analysis than OpenFOAM's native tools.

•	[**Porous Media Simulation with Energy Consideration**]: Facilitates cutting-edge porous media simulations, incorporating functionalities that simplify the process, including the explicit treatment of the energy equation.

•	[**Built-In Energy Equation Solutions**]: Unlike many OpenFOAM solvers that don't natively address the energy equation, SolarEnergy4Foam has this capability built into its solvers, obviating the need for extra coding.

•	[**Inclusion of Radiation Elements**]: Enhances solvers with radiation elements, enabling simulations that involve thermal radiation components.

•	[**Solvers for Incompressible Conjugate Heat Transfer**]: Incorporates solvers that handle incompressible conjugate heat transfer under both steady and transient conditions, as well as under forced and free convection scenarios.

With these attributes, SolarEnergy4Foam positions itself as an all-encompassing, intuitive toolbox that streamlines the CFD simulation workflow for solar energy applications, delivering features that boost both the precision and efficiency of research in this crucial sector.


## Getting Started
To utilize the SolarEnergy4Foam toolbox, follow these steps:

1. Install [**FoamExtend 4.1**].
2. Install the [**coolProp**] repository.
3. Clone the [**SE4Foam**] repository from GitHub.
4. Navigate to the SE4Foam directory and execute the `Allwmake` script to compile the toolbox.

## Contribution
We welcome contributions from the community. If you're interested in improving SE4Foam or have suggestions, please reach out to the authors or submit a pull request.


### How can I cite the current miniFoam and find more information?

You can see the governing equations and boundary conditions for each case in tutorial and more case in following resources:

[1] Vahedi, B., Golab, E., Sadr, A.N. and Vafai, K., 2022. Thermal, thermodynamic and exergoeconomic investigation of a parabolic trough collector utilizing nanofluids. Applied Thermal Engineering, 206, p.118117.

[2] Golab, E., Vahedi, B., Jain, A., Taylor, R.A. and Vafai, K., 2023. Laminar forced convection in a tube with a nano-encapsulated phase change materials: Minimizing exergy losses and maximizing the heat transfer rate. Journal of Energy Storage, 65, p.107233.

## License
SE4Foam is released under an open-source license. For more details, see the LICENSE file in the repository.

## Contact
For any inquiries or support, please contact the authors via the provided email addresses.

---

Embrace the future of solar energy with SolarEnergy4Foam, your partner in sustainable energy solutions.


