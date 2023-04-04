
REM arguments are in the following order:  SD_DELTA_CYCLE, P_DEATH_T_ZERO, EC50, LINEAGE_WEIGHT, ALIEN_WEIGHT, RSEED
echo start simulation without competition
SimuCellWinx64.exe 250 0.0025 60 1 1 11111

echo start simulation with competition
SimuCellWinx64.exe 250 0.0025 60 0.1 3 11111