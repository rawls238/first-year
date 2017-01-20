
%% Load data for plotting
two_a_r_6 = [0.00108225,0.038961,0.243506,0.4329,0.243506,0.038961,0.00108225]; %r = 6 stationary distribution
two_a_r_7 = [0.000291375,0.0142774,0.128497,0.356935,0.356935,0.128497,0.0142774,0.000291375]; %r = 7 stationary distribution

simb_stock_1 = [100.0,106.0,112.36,119.102,126.248,133.823,141.852,150.363,159.385,168.948,179.085,189.83,201.22,213.293,226.09,239.656,254.035,269.277,285.434,302.56,320.714,339.956,360.354,381.975,404.893,429.187,454.938,482.235,511.169,541.839,574.349,608.81,645.339,684.059,725.103,768.609,814.725,863.609,915.425,970.351,1028.57,1090.29,1155.7,1225.05,1298.55,1376.46,1459.05,1546.59,1639.39,1737.75,1842.02,1952.54,2069.69,2193.87,2325.5,1860.4,1972.03,2090.35,2215.77,2348.71,2489.64,2639.02,2797.36,2965.2,3143.11,3331.7,3531.6,2825.28,2260.22,1808.18,1446.54,1157.23,1226.67,1300.27,1378.28,1460.98,1548.64,1641.56,1740.05,1844.46,1475.56,1180.45,1251.28,1001.02,800.818,848.867,899.799,953.787,1011.01,1071.68,1135.98,1204.13,1276.38,1352.97,1434.14,1520.19,1611.4,1708.09,1366.47,1448.46,1535.37];
simb_stock_2 = [100.0,106.0,112.36,119.102,126.248,133.823,141.852,150.363,159.385,168.948,179.085,189.83,201.22,213.293,226.09,239.656,254.035,269.277,285.434,302.56,320.714,339.956,360.354,381.975,404.893,429.187,454.938,482.235,511.169,541.839,574.349,608.81,645.339,684.059,725.103,768.609,814.725,863.609,915.425,970.351,1028.57,1090.29,1155.7,924.563,980.036,1038.84,1101.17,1167.24,1237.27,1311.51,1390.2,1473.61,1562.03,1655.75,1755.1,1860.4,1972.03,2090.35,2215.77,2348.71,2489.64,2639.02,2797.36,2965.2,3143.11,3331.7,3531.6,3743.49,3968.1,4206.19,4458.56,4726.07,5009.64,5310.22,5628.83,5966.56,6324.55,6704.03,7106.27,7532.64,7984.6,8463.68,8971.5,9509.79,10080.4,10685.2,11326.3,12005.9,12726.2,13489.8,14299.2,15157.2,16066.6,17030.6,18052.4,19135.6,20283.7,21500.7,22790.8,24158.2,25607.7];
simb_stock_3 = [100.0,106.0,112.36,119.102,126.248,133.823,141.852,150.363,159.385,168.948,179.085,189.83,201.22,213.293,226.09,239.656,254.035,269.277,285.434,302.56,320.714,339.956,360.354,381.975,404.893,429.187,454.938,482.235,511.169,541.839,574.349,608.81,645.339,684.059,725.103,768.609,814.725,863.609,915.425,970.351,1028.57,1090.29,1155.7,1225.05,1298.55,1376.46,1459.05,1546.59,1639.39,1737.75,1842.02,1952.54,2069.69,2193.87,2325.5,2465.03,2612.93,2769.71,2935.89,3112.05,3298.77,3496.7,2797.36,2965.2,3143.11,3331.7,3531.6,3743.49,3968.1,3174.48,3364.95,3566.85,3780.86,3024.69,3206.17,3398.54,3602.45,3818.6,4047.71,4290.58,4548.01,4820.89,5110.15,5416.75,5741.76,6086.27,6451.44,6838.53,7248.84,7683.77,8144.8,8633.48,9151.49,9700.58,10282.6,10899.6,11553.5,12246.8,12981.6,13760.5,14586.1];
simb_stock_4 = [100.0,106.0,112.36,119.102,126.248,133.823,107.058,113.482,120.29,127.508,135.158,143.268,151.864,160.976,170.634,180.872,191.725,203.228,215.422,228.347,242.048,256.571,271.965,288.283,305.58,323.915,343.35,363.951,385.788,408.935,433.471,459.479,487.048,516.271,547.247,580.082,614.887,651.78,690.887,732.34,776.281,822.857,872.229,924.563,980.036,1038.84,1101.17,1167.24,1237.27,1311.51,1390.2,1473.61,1562.03,1655.75,1755.1,1860.4,1972.03,2090.35,2215.77,1772.61,1878.97,1991.71,2111.21,2237.88,2372.16,2514.49,2665.36,2825.28,2994.79,3174.48,3364.95,3566.85,3780.86,4007.71,4248.17,4503.06,4773.25,5059.64,5363.22,5685.01,6026.12,6387.68,6770.94,7177.2,7607.83,8064.3,8548.16,9061.05,9604.71,10181.0,10791.9,11439.4,12125.7,12853.3,13624.5,14441.9,15308.5,16227.0,17200.6,18232.6,19326.6];

simb_liq_1 = [1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,1,1,1,1,1,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1];
simb_liq_2 = [1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1];
simb_liq_3 = [1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1];
simb_liq_4 = [1,1,1,1,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1];

simc_stock_1 = [100.0,106.0,112.36,119.102,126.248,133.823,141.852,150.363,159.385,168.948,179.085,189.83,201.22,213.293,226.09,239.656,254.035,269.277,285.434,302.56,320.714,339.956,360.354,381.975,404.893,429.187,454.938,482.235,511.169,541.839,574.349,608.81,645.339,684.059,725.103,768.609,814.725,651.78,690.887,732.34,776.281,822.857,658.286,526.629,558.226,591.72,627.223,664.857,704.748,563.798,597.626,633.484,671.493,537.194,569.426,603.592,482.873,386.299,409.477,327.581,262.065,209.652,222.231,235.565,249.699,264.681,280.562,297.395,315.239,252.191,201.753,161.402,171.087,181.352,192.233,153.786,123.029,130.411,138.235,110.588,117.224,124.257,99.4056,79.5245,63.6196,50.8957,53.9494,57.1864,45.7491,36.5993,29.2794,31.0362,32.8984,26.3187,27.8978,22.3183,23.6574,25.0768,20.0614,16.0491,17.0121];
simc_stock_2 = [100.0,106.0,112.36,89.888,95.2813,100.998,107.058,113.482,120.29,127.508,135.158,143.268,151.864,160.976,170.634,180.872,191.725,203.228,215.422,228.347,242.048,256.571,271.965,288.283,305.58,323.915,343.35,363.951,385.788,408.935,433.471,459.479,487.048,389.638,413.017,437.798,464.066,491.91,521.424,552.71,585.872,621.024,658.286,697.783,739.65,784.029,831.071,880.935,933.791,989.819,1049.21,1112.16,1178.89,1249.62,1324.6,1404.08,1488.32,1577.62,1262.1,1337.82,1418.09,1503.18,1202.54,962.033,1019.76,1080.94,1145.8,1214.54,1287.42,1364.66,1446.54,1533.33,1226.67,981.334,1040.21,832.172,882.102,705.681,564.545,451.636,478.734,507.458,405.967,324.773,259.819,275.408,220.326,176.261,141.009,112.807,90.2456,95.6604,101.4,81.12,85.9872,68.7898,72.9171,77.2922,81.9297,86.8455,92.0562];
simc_stock_3 = [100.0,106.0,112.36,119.102,126.248,133.823,141.852,150.363,120.29,127.508,135.158,108.127,114.614,91.6914,73.3531,58.6825,46.946,49.7628,39.8102,42.1988,44.7307,47.4146,37.9317,40.2076,32.1661,25.7328,27.2768,28.9134,23.1307,24.5186,25.9897,20.7918,22.0393,23.3616,18.6893,14.9514,15.8485,16.7994,13.4395,14.2459,15.1007,16.0067,16.9671,13.5737,14.3881,15.2514,12.2011,9.7609,10.3466,10.9673,8.77388,7.0191,7.44025,5.9522,6.30933,6.68789,7.08916,5.67133,4.53706,4.80929,5.09785,5.40372,4.32297,3.45838,2.7667,2.93271,2.34616,1.87693,1.50155,1.59164,1.27331,1.34971,1.07977,0.863814,0.915642,0.970581,1.02882,0.823053,0.872436,0.697949,0.739826,0.59186,0.473488,0.501898,0.532012,0.425609,0.340487,0.27239,0.217912,0.17433,0.184789,0.195877,0.207629,0.220087,0.17607,0.186634,0.149307,0.119446,0.0955565,0.0764452,0.0611562];
simc_stock_4 = [100.0,106.0,112.36,119.102,126.248,100.998,107.058,113.482,120.29,127.508,135.158,143.268,151.864,160.976,170.634,180.872,191.725,203.228,215.422,228.347,242.048,256.571,271.965,288.283,305.58,323.915,343.35,363.951,385.788,408.935,433.471,459.479,487.048,516.271,547.247,437.798,464.066,491.91,393.528,417.139,442.168,468.698,496.82,526.629,558.226,591.72,627.223,664.857,704.748,747.033,791.855,839.366,889.728,943.112,754.49,799.759,847.744,898.609,952.526,1009.68,1070.26,1134.47,1202.54,1274.69,1019.76,1080.94,864.753,916.638,733.31,777.309,823.947,873.384,698.707,558.966,592.504,628.054,502.443,401.955,426.072,451.636,478.734,382.987,405.967,324.773,259.819,207.855,166.284,133.027,141.009,149.469,158.437,126.75,101.4,81.12,85.9872,91.1464,72.9171,77.2922,61.8337,65.5438,52.435];

simc_liq_1 = [1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0];
simc_liq_2 = [1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0];
simc_liq_3 = [1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0];
simc_liq_4 = [1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0];

%% question 2
figure()
subplot(1,2,1)
bar(two_a_r_6);
title('r = 6')
subplot(1,2,2)
bar(two_a_r_7)
title('r = 7');
%% Quesiton 4, part b
figure()
subplot(2,2,1);
hold on; box on;
plot(1:1:100, simb_stock_1(1:100));
yyaxis right;
plot(1:1:100, simb_liq_1(1:100));
ylim([-1 2]);
l=legend('Stock Price','Liquidity', 'Location','NorthWest');
legend('boxoff')

subplot(2,2,2);
hold on; box on;
plot(1:1:100, simb_stock_2(1:100));
yyaxis right;
plot(1:1:100, simb_liq_2(1:100));
ylim([-1 2]);
l=legend('Stock Price','Liquidity', 'Location','NorthWest');
legend('boxoff')

subplot(2,2,3);
hold on; box on;
plot(1:1:100, simb_stock_3(1:100));
yyaxis right;
plot(1:1:100, simb_liq_3(1:100));
ylim([-1 2]);
l=legend('Stock Price','Liquidity', 'Location','NorthWest');
legend('boxoff')

subplot(2,2,4);
hold on; box on;
plot(1:1:100, simb_stock_4(1:100));
yyaxis right;
plot(1:1:100, simb_liq_4(1:100));
ylim([-1 2]);
l=legend('Stock Price','Liquidity', 'Location','NorthWest');
legend('boxoff')

% part c
figure()
subplot(2,2,1);
hold on; box on;
plot(1:1:100, simc_stock_1(1:100));
yyaxis right;
plot(1:1:100, simc_liq_1(1:100));
ylim([-1 2]);
l=legend('Stock Price','Liquidity', 'Location','NorthWest');
legend('boxoff')

subplot(2,2,2);
hold on; box on;
plot(1:1:100, simc_stock_2(1:100));
yyaxis right;
plot(1:1:100, simc_liq_2(1:100));
ylim([-1 2]);
l=legend('Stock Price','Liquidity', 'Location','NorthWest');
legend('boxoff')

subplot(2,2,3);
hold on; box on;
plot(1:1:100, simc_stock_3(1:100));
yyaxis right;
plot(1:1:100, simc_liq_3(1:100));
ylim([-1 2]);
l=legend('Stock Price','Liquidity', 'Location','NorthWest');
legend('boxoff')

subplot(2,2,4);
hold on; box on;
plot(1:1:100, simc_stock_4(1:100));
yyaxis right;
plot(1:1:100, simc_liq_4(1:100));
ylim([-1 2]);
l=legend('Stock Price','Liquidity', 'Location','NorthWest');
legend('boxoff')