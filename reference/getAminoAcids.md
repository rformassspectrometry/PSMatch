# Amino acids

Returns a `data.frame` of amino acid properties: `AA`, `ResidueMass`,
`Abbrev3`, `ImmoniumIonMass`, `Name`, `Hydrophobicity`,
`Hydrophilicity`, `SideChainMass`, `pK1`, `pK2` and `pI`.

## Usage

``` r
getAminoAcids()
```

## Value

`data.frame`

## Author

Laurent Gatto

## Examples

``` r

getAminoAcids()
#>     AA ResidueMass Abbrev3 ImmoniumIonMass                Name Hydrophobicity
#> 1  peg    44.00000    <NA>              NA Polyethylene glycol             NA
#> 2    A    71.03711     Ala        44.05003             Alanine           0.62
#> 3    R   156.10111     Arg       129.11400            Arginine          -2.53
#> 4    N   114.04293     Asn        87.05584          Asparagine          -0.78
#> 5    D   115.02694     Asp        88.03986       Aspartic acid          -0.90
#> 6    C   103.00919     Cys        76.02210            Cysteine           0.29
#> 7    E   129.04259     Glu       102.05550       Glutamic acid          -0.74
#> 8    Q   128.05858     Gln       101.07150           Glutamine          -0.85
#> 9    G    57.02146     Gly        30.03438             Glycine           0.48
#> 10   H   137.05891     His       110.07180           Histidine          -0.40
#> 11   I   113.08406     Ile        86.09698          Isoleucine           1.38
#> 12   L   113.08406     Leu        86.09698             Leucine           1.06
#> 13   K   128.09496     Lys       101.10790              Lysine          -1.50
#> 14   M   131.04049     Met       104.05340          Methionine           0.64
#> 15   F   147.06841     Phe       120.08130       Phenylalanine           1.19
#> 16   P    97.05276     Pro        70.06568             Proline           0.12
#> 17   S    87.03203     Ser        60.04494              Serine          -0.18
#> 18   T   101.04768     Thr        74.06059           Threonine          -0.05
#> 19   W   186.07931     Trp       159.09220          Tryptophan           0.81
#> 20   Y   163.06333     Tyr       136.07620            Tyrosine           0.26
#> 21   V    99.06841     Val        72.08133              Valine           1.08
#> 22   U   150.95363     Sec       122.95860      Selenocysteine             NA
#> 23   O   237.14772     Pyl       210.16060         Pyrrolysine             NA
#>    Hydrophilicity SideChainMass  pK1   pK2    pI
#> 1              NA            NA   NA    NA    NA
#> 2            -0.5            15 2.35  9.87  6.11
#> 3             3.0           101 2.18  9.09 10.76
#> 4             0.2            58 2.18  9.09 10.76
#> 5             3.0            59 1.88  9.60  2.98
#> 6            -1.0            47 1.71 10.78  5.02
#> 7             3.0            73 2.19  9.67  3.08
#> 8             0.2            72 2.17  9.13  5.65
#> 9             0.0             1 2.34  9.60  6.06
#> 10           -0.5            82 1.78  8.97  7.64
#> 11           -1.8            57 2.32  9.76  6.04
#> 12           -1.8            57 2.36  9.60  6.04
#> 13            3.0            73 2.20  8.90  9.47
#> 14           -1.3            75 2.28  9.21  5.74
#> 15           -2.5            91 2.58  9.24  5.91
#> 16            0.0            42 1.99 10.60  6.30
#> 17            0.3            31 2.21  9.15  5.68
#> 18           -0.4            45 2.15  9.12  5.60
#> 19           -3.4           130 2.38  9.39  5.88
#> 20           -2.3           107 2.20  9.11  5.63
#> 21           -1.5            43 2.29  9.74  6.02
#> 22             NA            94 2.19  9.21  7.25
#> 23             NA           198 2.20  9.00  9.75
```
