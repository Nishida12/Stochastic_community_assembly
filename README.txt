# N_Species_M_Resources_Figure2A.cpp and N_Species_M_Resources_Figure2B.cpp

## description
(N_Species_M_Resources_Figure2A.cpp)
All resources have some microbial species specialized for each resource.
Each microorganism takes up resources and secretes metabolic by-products.
(N_Species_M_Resources_Figure2B.cpp)
Microbial communities adapted to high/low resource concentrations have been supplied a resource of high/low concentration.
The communities adapted to a high resource concentration consist of a large number of species adapted to high resource concentration and a small number of species adapted to low resource concentration.
The communities adapted to a low resource concentration have the opposite composition.

## preparation
Make folder named "Result", folder named "TimeLapse" into "Result" folder, and 100 folders named "0", "1", ... , "99" into "TimeLapse" folder.

## parameter setting in N_Species_M_Resources_Figure2B.cpp
For supplying a high resource concentration, supply resource concentration K becomes 100.
For supplying a low resource concentration, supply resource concentration K becomes 1.
For communities adapted to a high resource concentration,
C_maj = 4.0 * K / (K_h + K) * (double)K_h / (K_h + K)
C_min = 4.0 * K / (K_l + K) * (double)K_l / (K_l + K)
For communities adapted to a low resource concentration,
C_maj = 4.0 * K / (K_l + K) * (double)K_l / (K_l + K)
C_min = 4.0 * K / (K_h + K) * (double)K_h / (K_h + K)
