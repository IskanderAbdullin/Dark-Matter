# Input for the independent kinetic benchmarks

Use the [Saikawa–Shirai 2020 central equation-of-state table](https://member.ipmu.jp/satoshi.shirai/DM2020/eos2020.dat). This external input is not committed to the repository. From the repository root, run:

```sh
cd kinetics
curl -fsSL https://member.ipmu.jp/satoshi.shirai/DM2020/eos2020.dat -o EOS.full.dat
LC_ALL=C awk 'NF == 3 && $1 >= 0.009 && $1 <= 1.1 {print $1, $2, $3}' EOS.full.dat > EOS.dat
printf '%s\n' 'e19bdd99a4891a9b1e14c51f4831d32abf9e7a8c4cd144fb227e29e93cd6f4bd  EOS.dat' | sha256sum -c -
```

The benchmark file contains all 98 table rows in the requested range
`0.009 <= T [GeV] <= 1.1`, with original numeric tokens separated by
single spaces and one final newline. The temperature endpoints in the
discrete file are 0.00925767 and 1.05163 GeV. The SHA-256 check pins this
filtered input without redistributing the data. Python's `StandardModel`
class applies the same temperature filter to a supplied full table.

For the late tabulated radiation clock, use the original temperature
coverage in `EOS.full.dat`, not the 98-row warm input. The verified
canonical extraction for the late runs used every 3-column row with
`0 <= T [GeV] <= 1.1` (zero-temperature row retained for verification;
not interpolated). Optional check:

```sh
LC_ALL=C awk 'NF == 3 && $1 >= 0 && $1 <= 1.1 {print $1, $2, $3}' EOS.full.dat > EOS.late.dat
printf '%s\n' 'bc95d6e13617c5ae8078c6f7bf2f8323beee2b8f4aee56082c7635fa038a238f  EOS.late.dat' | sha256sum -c -
```
