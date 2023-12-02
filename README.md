# Nechako-sturgeon-spatial-abundance
Brownie multi-state estimation model for abundance

R code simulates spatial movement of fish through a river network. Simulated fish are captured and have tags applied. A subset of fish have telemetry tags which are used to determine movement. 

R code then calls TMB code which can estimate movement, abundance, and mortality parameters. 

This code package can also directly estimate parameters from real data by setting sim==FALSE and providing files "Nechako tagging.csv" and "Nechako telemetry.csv"
