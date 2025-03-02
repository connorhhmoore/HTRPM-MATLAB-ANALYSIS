Script to analyze sphere packing in cylindrical cores of pebble-bed reactors. Example output using HTR-PM core geometry and 6cm pebbbles below:

```
 Calculated values using a core D/H (3.00 [m])/(11.00 [m]).
 Fuel pebble diameter 0.06 [m].
 Calculation time of 4.813291e+02 seconds with ext = 4

 ======== Basic stacking (box/cube stacking) pattern ========
  -> Total number of contained pellets: 344955 [n]
  -> Total surface area of all pellets: 3901.35 [m²]
  -> Calculated surface heat flux of the core: 64080.40 [W/m²]
  -> Average number of pellets per cubic meter: 4436.47 [n/m³]
  -> Average number of pellets per sheet: 1885.00 [n/sheet]
  -> Volume packing efficiency: 50.18 [%]

 ======== BCC-style stacking (ABAB sheet-based stacking) pattern ========
  -> Total number of contained pellets: 421232 [n]
  -> Total surface area of all pellets: 4764.02 [m²]
  -> Calculated surface heat flux of the core: 52476.67 [W/m²]
  -> Average number of pellets per cubic meter: 5417.47 [n/m³]
  -> Average number of pellets per sheet: 1880.50 [n/sheet]
  -> Volume packing efficiency: 61.27 [%]

 ======== HCP-style stacking (ABCABC sheet-based stacking) pattern ========
  -> Total number of contained pellets: 487501 [n]
  -> Total surface area of all pellets: 5513.51 [m²]
  -> Calculated surface heat flux of the core: 45343.20 [W/m²]
  -> Average number of pellets per cubic meter: 6269.75 [n/m³]
  -> Average number of pellets per sheet: 2176.34 [n/sheet]
  -> Volume packing efficiency: 70.91 [%]

 ======== Real analysis (no stacking calculation) ========
  -> Total claimed number of pellets: 425000 [n]
  -> Volume packing efficiency: 61.82 [%]
  -> Total surface area of all pellets: 4806.64 [m²]
  -> Calculated surface heat flux of the core: 52011.42 [W/m²]
  -> ΔT to center of fuel pellet: 0.002 [°C]

 ======== Flow cross-section calculations for BCC ========
  -> Center-to-center sphere spacing: 0.049 [m]
  -> XS flow area for entire reactor core 7.069 [m²]
  -> Maximum XS area for one pebble: 0.003 [m²]
  -> Maximum XS area for reactor core: 3.872 [m²]
  -> Minimum XS area for reactor core: 1.739 [m²]
  -> Average XS flow area for reactor core 2.717 [m²]

 (Sphere packing analysis, Connor Moore, 2024-2025)
      <connor.moore@ontariotechu.net>
 (Calculations performed 02-Mar-2025 14:34:25)
```
