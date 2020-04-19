# vcf-to-wegene

### Uasge
```
usage: vcf_to_wegene.py [-h] -i INPUT -b BLANK [-o OUTPUT]

Convert allsite VCF to wegene/23andMe raw data format

optional arguments:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        allsite.vcf.gz with index [requested]
  -b BLANK, --blank BLANK
                        wegene/23andMe format file without GT [requested]
  -o OUTPUT, --output OUTPUT
                        output wegene/23andMe format file
```

### 示例:
```
python vcf_to_wegene.py -i input.allsite.vcf.gz -b wegene_blank.txt.gz -o out.wegene.txt
```
