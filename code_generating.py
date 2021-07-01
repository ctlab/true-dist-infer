vs = range(50, 2001, 10)

for v in vs:
    print(f'cp -r ucsc_7_working ucsc_7_diff_res/{v}k')
    print(f'sed -i -e \'s/50000/{v * 1000}/g\' ucsc_7_diff_res/{v}k/config.file')
