## pQuant standard file format

We have defined and standard format for all packages output to provide the final outputs including Mstats, Proteus or Triqler.

First a configuration file is needed to describe `contrasts`, samples, fractions, etc. The configuration file:

```json
{
  "configuration": {
    "analysis": {
      "assay_groups": {
      "assay_group": [
          {
            "assay": [
              "4.AD_B",
              "4.AD_A",
              "4.AD_C"
            ]
          },
          {
            "assay": [
              "3.AD_B",
              "3.AD_C",
              "3.AD_A"
            ]
          },
          {
            "assay": [
              "2.AD_B",
              "2.AD_A",
              "2.AD_C"
            ]
          },
          {
            "assay": [
              "1.AD_B",
              "1.AD_A",
              "1.AD_C"
            ]
          }
        ]
      },
      "contrasts": {
        "contrast": [
          {
            "" A": "'Contrast 1' vs 'Contrast 1",
            "reference_assay_group": "g2",
            "test_assay_group": "g1"
          },
          {
            "name": "'Alzheimer's disease' vs 'normal' in 'sexagenarian age bracket'",
            "reference_assay_group": "g4",
            "test_assay_group": "g3"
          }
        ]
      }
    }
  }
}
```
