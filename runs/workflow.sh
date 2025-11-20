python scripts\01_afnd_converter.py 

python scripts\02_select_hla_alleles.py ^
    --freq-table C:\Users\Kacper\Desktop\immuno\outputs\01_afnd\hla_freq.tsv ^
    --population "Cambodia" "Indonesia Java Western" "Indonesia Java pop 2" "Indonesia Sundanese and Javanese" "Indonesia Java" "Indonesia Moluccan Islands" "Indonesia Nusa Tenggara Islands" "Indonesia Java Yogyakarta Region" "Indonesia Singaporean" "Malaysia Champa" "Malaysia Kelantan" "Malaysia Mandailing" "Malaysia Patani" "Malaysia Peninsular Chinese" "Malaysia Peninsular Indian" "Malaysia Peninsular Malay" "Malaysia Jelebu Temuan" "Malaysia Kedah Baling Kensiu" "Malaysia Kedah Kensiu" "Malaysia Pahang Semai" "Malaysia Perak Grik Jehai" "Malaysia Sarawak Bau Bidayuh" "Malaysia Sabah Dusun" "Malaysia Sabah Kadazan" "Malaysia Sarawak Bidayuh" "Malaysia Sarawak Iban" "Malaysia" "Philippines Ivatan" "Philippines" "Singapore Chinese" "Singapore Chinese Han" "Singapore Javaneses" "Singapore Riau Malay" "Singapore SGVP Chinese CHS" "Singapore SGVP Malay MAS" "Singapore SGVP. Indian INS" "Singapore Thai" "Thailand" "Thailand Northeast" "Thailand Northeast pop 2" "Thailand Bangkok" "Thailand East Chanthaburi Khmer" "Thailand Kamphaeng Phet" "Thailand North Dai Lue" "Thailand pop 3" "Thailand pop 2" "Thailand Bangkok pop 2" "Thailand Bangkok pop 3" "Vietnam Hanoi Kinh pop 2" "Vietnam Hanoi Kinh" "Vietnam HoaBinh Muong" "Vietnam Kinh" "Vietnam Kinh DQB1" ^
    --netmhcpan-alleles C:\Users\Kacper\Desktop\immuno\data\MHC_allele_names.txt ^
    --netmhciipan-alleles C:\Users\Kacper\Desktop\immuno\data\alleles_name.txt ^
    --max-classI 20 ^
    --max-classII 15 ^
    --coverage 1.0
