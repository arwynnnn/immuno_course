## venv
### linux
```bash
python3 -m venv .venv
source .venv/bin/activate
pip install pandas requests tqdm
```

## windows
```bash
python3 -m venv .venv
source .venv\Scripts\activate
pip install pandas requests tqdm
```

# Population selection
As per [https://doi.org/10.1002/cncr.28312](https://doi.org/10.1002/cncr.28312) and [https://pmc.ncbi.nlm.nih.gov/articles/PMC5411279](https://pmc.ncbi.nlm.nih.gov/articles/PMC5411279):
> There is a tremendous discrepancy in the worldwide incidence of ICC, with significantly higher rates of ICC seen in Eastern Asia when compared to Western countries

We will choose those populations for selection HLA allels, aiming for the broadest coverege for those:

```
'Cambodia' 'Indonesia Java Western' 'Indonesia Java pop 2' 'Indonesia Sundanese and Javanese' 'Indonesia Java' 'Indonesia Moluccan Islands' 'Indonesia Nusa Tenggara Islands' 'Indonesia Java Yogyakarta Region' 'Indonesia Singaporean' 'Malaysia Champa' 'Malaysia Kelantan' 'Malaysia Mandailing' 'Malaysia Patani' 'Malaysia Peninsular Chinese' 'Malaysia Peninsular Indian' 'Malaysia Peninsular Malay' 'Malaysia Jelebu Temuan' 'Malaysia Kedah Baling Kensiu' 'Malaysia Kedah Kensiu' 'Malaysia Pahang Semai' 'Malaysia Perak Grik Jehai' 'Malaysia Sarawak Bau Bidayuh' 'Malaysia Sabah Dusun' 'Malaysia Sabah Kadazan' 'Malaysia Sarawak Bidayuh' 'Malaysia Sarawak Iban' 'Malaysia' 'Philippines Ivatan' 'Philippines' 'Singapore Chinese' 'Singapore Chinese Han' 'Singapore Javaneses' 'Singapore Riau Malay' 'Singapore SGVP Chinese CHS' 'Singapore SGVP Malay MAS' 'Singapore SGVP. Indian INS' 'Singapore Thai' 'Thailand' 'Thailand Northeast' 'Thailand Northeast pop 2' 'Thailand Bangkok' 'Thailand East Chanthaburi Khmer' 'Thailand Kamphaeng Phet' 'Thailand North Dai Lue' 'Thailand pop 3' 'Thailand pop 2' 'Thailand Bangkok pop 2' 'Thailand Bangkok pop 3' 'Vietnam Hanoi Kinh pop 2' 'Vietnam Hanoi Kinh' 'Vietnam HoaBinh Muong' 'Vietnam Kinh' 'Vietnam Kinh DQB1'
```

## output

We take the precomputed up-to-date frequency table generated using this [repo](https://github.com/slowkow/allelefrequencies/tree/main).

Moreover, we limit allels to those available to use at the DTU Health:
- [NetMHCpan-4.1](https://services.healthtech.dtu.dk/services/NetMHCpan-4.1/MHC_allele_names.txt)
- [NetMHCIIpan-4.3](https://services.healthtech.dtu.dk/services/NetMHCIIpan-4.3/alleles_name.txt)

```
Selected Class I alleles (n=13):
   locus       allele  frequency netmhcpan_allele
0      A  HLA-A*24:07     0.5000       HLA-A24:07
1      A  HLA-A*02:01     0.4050       HLA-A02:01
2      A  HLA-A*11:01     0.3400       HLA-A11:01
3      B  HLA-B*18:01     0.2800       HLA-B18:01
4      B  HLA-B*15:13     0.2780       HLA-B15:13
5      B  HLA-B*15:02     0.2200       HLA-B15:02
6      B  HLA-B*40:01     0.2040       HLA-B40:01
7      B  HLA-B*13:01     0.2000       HLA-B13:01
8      C  HLA-C*08:01     0.3100       HLA-C08:01
9      C  HLA-C*07:02     0.2230       HLA-C07:02
10     C  HLA-C*07:01     0.1840       HLA-C07:01
11     C  HLA-C*01:02     0.1753       HLA-C01:02
12     C  HLA-C*06:02     0.1730       HLA-C06:02

Selected Class II alleles (n=8):
  locus          allele  frequency netmhciipan_allele
0  DRB1  HLA-DRB1*07:01      0.903         DRB1*07:01
1  DRB1  HLA-DRB1*12:02      0.534         DRB1*12:02
2  DQB1  HLA-DQB1*03:01      0.494         DQB1*03:01
3  DQB1  HLA-DQB1*05:02      0.481         DQB1*05:02
4  DQB1  HLA-DQB1*03:03      0.313         DQB1*03:03
5  DPB1  HLA-DPB1*05:01      0.500         DPB1*05:01
6  DPB1  HLA-DPB1*04:01      0.428         DPB1*04:01
7  DPB1  HLA-DPB1*01:01      0.337         DPB1*01:01
```