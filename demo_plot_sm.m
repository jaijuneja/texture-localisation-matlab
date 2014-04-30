num_maps = ...
[1
4
8
16
24
36
54
80
111
170
356
548
646];

branching = ...
[17.3423
9.2266
6.4418
4.9373
4.4698
4.0239
3.0772
2.8757
2.5059
2.03
1.6736
1.4725
1.3194];

overlap = ...
[0.93858
0.93851
0.9412
0.94772
0.94757
0.93688
0.93438
0.93952
0.92863
0.8905
0.8713
0.94858
0.30228];
overlap = overlap * 100;

time = ...
[1
0.878402339
0.724067874
0.492037538
0.241055626
0.227434886
0.07386893
0.051388986
0.034051985
0.016029207
0.008279655
0.005015108
0.002564258];
time = time * 100;
% 
figure;
[haxes,hline1,hline2] = plotyy(num_maps,branching,num_maps,overlap);
ylabel(haxes(1),'Avg. Branching Factor per Submap') % label left y-axis
ylabel(haxes(2),'Localisation Accuracy (% of Baseline)') % label right y-axis
xlabel(haxes(2),'Number of Submaps') % label x-axis
set(hline1, 'LineWidth', 2)
set(hline2, 'LineWidth', 2)
set(haxes(2), 'ylim', [0 100], 'FontSize', 14)
set(haxes(1), 'ylim', [0 20], 'FontSize', 14)
set(haxes(1), 'ytick', [0 4 8 12 16 20])

figure;
plot(branching, time, 'LineWidth', 2)
xlabel('Avg. Branching Factor per Submap')
ylabel('Computation Time (% of Baseline)')
set(gca, 'FontSize', 14)