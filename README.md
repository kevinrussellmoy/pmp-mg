# pmp-mg
Implementation of Pontryagin's Minimum Principle for microgrid energy storage control.

![mg_diagram2](https://user-images.githubusercontent.com/68761709/206780848-cc2e6ca9-be52-448f-8405-763cf7723429.png)

![fuelcons_degrad_output_phase1_SOCs](https://user-images.githubusercontent.com/68761709/206780982-f9adfe9b-7ce5-4a0e-9c9a-23e88cd0a628.png)

This work resulted in a conference paper publication at the American Control Conference (ACC) 2022 in Atlanta, Georgia.

Abstract: Microgrids are energy systems that are able to supply power reliably in the face of instability on the main electric grid, increasingly driven by the effects of anthropogenic climate change. Microgrids are powered by diesel generators, energy storage, and renewable energy resources such as photovoltaics, to supply power to loads. Lithium-ion batteries (LIBs) are currently the dominant grid-scale energy storage technology and leading candidate for deployment in microgrids. An optimal control problem can be formulated regarding the optimal energy management of the LIB and other microgrid components, with the goal of minimizing the fuel consumption of the diesel engine. In this paper, Pontryagin's Minimum Principle (PMP) is used to solve the optimal energy management problem where the LIB is modeled through an equivalent circuit model. A semi-empirical model is used to assess the degradation of the LIB under the resulting optimal control. PMP is applied to a variety of initial and final outage conditions taken from real-world scenarios, resulting in different impacts on LIB degradation and fuel consumption.

Full conference paper available here: https://pangea.stanford.edu/ERE/pdf/OnoriPDF/Conferences/91.pdf
