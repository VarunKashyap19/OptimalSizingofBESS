#Hierarchical Optimal Allocation of Battery Energy Storage Systems for Multiple Services in Distribution Systems

#🌟 Overview
This repository contains the research paper "Hierarchical Optimal Allocation of Battery Energy Storage Systems for Multiple Services in Distribution Systems." It introduces a hierarchical optimization framework to strategically allocate Battery Energy Storage Systems (BESS) for enhanced voltage stability and reduced power losses in distribution networks.

#🎯 Why is this important?
With the increasing integration of Distributed Generation (DG) and the unpredictability of power injections, the grid faces voltage violation risks and power instability. This work presents a two-level hierarchical approach:

High-Level Allocation: Optimal placement of BESS to handle voltage deviations due to power uncertainties.
Low-Level Operation: Real-time control of BESS to support grid stability and minimize power losses.

#🛠️ Core Methodology
Electrical Constraints:

Power Balance
Voltage and Current Limits
Substation Power Limit
BESS Operation Constraints:

Power Limitation
Energy Balance over the Planning Period
State of Charge (SOC) Management
Convexification via Relaxation: Using a matrix variable 𝑊𝑡 to encode voltage and power flow relationships, the problem is solved through Second-Order Cone Programming (SOCP) for efficient optimization.

#⚡ Key Equation
The voltage deviation due to power disturbances is expressed as:

Δ𝑉 = 𝑅Δ𝑃 + 𝑋Δ𝑄
Where:
ΔP and ΔQ are active and reactive power injection disturbances.
R and X are matrices reflecting network characteristics.

#🚀 Applications
Optimal BESS placement for voltage regulation
Improved grid resilience with renewable energy integration
Reduced power losses and enhanced stability.

#📈 Results
Significant reduction in voltage deviations
Enhanced power system reliability
Cost-effective BESS allocation strategy
