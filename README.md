# CanDriS
A Statistical Framework for Posterior Profiling Cancer-Driving Sites based on Recurrent Somatic Mutations

Contact
---
	Gu Xun			xgu@iastate.edu
	Zhou Zhan		zhanzhou@zju.edu.cn
		
Description
---
	CanDriS: A Statistical Framework for Posterior Profiling Cancer-Driving Sites based on Recurrent Somatic Mutations
	Version: V1.1
	Date: 2019.07.16
	Copyright 2017-2019. All Rights Reserved.
	
Usage
---	
	Example：
	perl candris.pl -r <reference> -i <ICGC(.tsv)/TCGA(.maf)/preprocessed_data>  -o1 <genes_with_driver_mutation_sites> -o2 <driver_mutation_sites> 
	
	Paramaters:
		-r	input ensembl reference file 
		-i	input mutation file	
		-o1	output file 1, mutation gene list 
		-o2	output file 2, mutation site list
		-h	help
	
	Perl Modules:
		Getopt::Long
		Statistics::Distributions
	
Analysis steps
---	
	1	Data preparation
		1.1 Choose specific version of Ensembl reference file. #ICGC(75_GRCh37.p13)/TCGA(75_GRCh38.p13)/preprocessed data# The first two are processed and offered beforehand.
		1.2 Collect mutation data from ICGC(.tsv)/TCGA(.maf)/preprocessed data and save it as an mutation file. Its format should be as follows:
			gene_id	trancript_id	DNA_mutation	protein_mutation	mutation_count
			ENSG00000084734 ENST00000264717 1740C>A  F580L  1
			ENSG00000120669 ENST00000379881  389G>A  S130N  1
			ENSG00000011258 ENST00000415868  547A>G  T183A  1
			ENSG00000160716 ENST00000368476  758T>C  I253T  1
			ENSG00000141867 ENST00000263377 3139C>T R1047W  1
		
	2	CanDriS input
		2.1 Ensembl reference file  
		2.2 Mutation file 
		
	3	Perform CanDriS
		As shown by usage.
		
	4	Output files in your path 
		You will get two output files: gene_with_driver_mutation_sites and driver_mutation_sites.
		Gene_with_driver_mutation_sites:
			Gene_id	Transcript_id	Chr	Gene_name	Max_hit	Count	Protein_length	Distribution	mut_site	f0	mean	variance	m0	m1	eta	LRT_p-value	Q(z)	zc	number_of_predited_sites
			ENSG00000113805	ENST00000263665	CNTN3	3	3	140	1028	3:417;2:76;2:93;2:111;2:223;2:455;2:508;2:617;2:636;2:641;2:741;2:811;2:851;2:962;1:30;1:31;1:32;1:40;1:43;1:57;1:66;1:69;1:74;1:81;1:87;1:107;1:115;1:124;1:129;1:135;1:157;1:159;1:160;1:161;1:171;1:174;1:189;1:197;1:222;1:234;1:259;1:263;1:270;1:277;1:278;1:300;1:311;1:322;1:336;1:356;1:361;1:365;1:381;1:383;1:391;1:408;1:426;1:438;1:439;1:442;1:446;1:461;1:467;1:499;1:519;1:524;1:526;1:539;1:587;1:591;1:605;1:608;1:612;1:621;1:624;1:628;1:640;1:643;1:646;1:656;1:662;1:664;1:667;1:668;1:675;1:686;1:693;1:696;1:701;1:709;1:713;1:714;1:725;1:731;1:733;1:735;1:740;1:745;1:774;1:779;1:783;1:797;1:803;1:805;1:807;1:814;1:816;1:837;1:850;1:864;1:882;1:886;1:898;1:902;1:910;1:912;1:924;1:930;1:937;1:963;1:964;1:969;1:975;1:980;1:1027;	125	0.8784046692607	0.136186770428016	0.148913195852073	0.129647892598125	2.08245733169802	0.0033484464479568	0.17218	0.663901767239021:3:1;0.109510497064927:2:13;0.00759808051741862:1:111;0.00047642954522531:0:903;
			ENSG00000110237	ENST00000263674	ARHGEF17	11	3	113	2063	3:459;2:263;2:849;2:1888;2:1962;1:226;1:227;1:278;1:280;1:283;1:287;1:305;1:338;1:347;1:348;1:364;1:379;1:413;1:415;1:441;1:445;1:461;1:468;1:473;1:485;1:512;1:519;1:531;1:537;1:542;1:556;1:560;1:570;1:586;1:604;1:628;1:631;1:647;1:650;1:665;1:682;1:692;1:707;1:717;1:721;1:722;1:778;1:831;1:834;1:875;1:895;1:918;1:931;1:949;1:964;1:966;1:1007;1:1033;1:1042;1:1067;1:1158;1:1166;1:1204;1:1224;1:1228;1:1256;1:1257;1:1307;1:1308;1:1326;1:1359;1:1363;1:1364;1:1400;1:1418;1:1419;1:1432;1:1434;1:1443;1:1579;1:1599;1:1604;1:1616;1:1619;1:1631;1:1683;1:1699;1:1702;1:1715;1:1741;1:1752;1:1771;1:1793;1:1807;1:1837;1:1883;1:1892;1:1896;1:1945;1:1957;1:1971;1:1974;1:1985;1:1999;1:2017;1:2024;1:2037;	107	0.948133785748909	0.0547746000969462	0.0585889768133098	0.0532596624764891	2.57261867820018	0.000601318673123638	0.057062	0.845189363222128:3:1;0.101547941724322:2:4;0.00233444789816806:1:102;4.84396684991006e-005:0:1956;
			ENSG00000145850	ENST00000274532	TIMD4	5	4	76	378	4:67;3:128;2:8;2:29;2:44;2:51;2:147;2:202;2:282;2:301;1:2;1:4;1:5;1:13;1:19;1:31;1:33;1:34;1:36;1:47;1:57;1:70;1:73;1:79;1:92;1:94;1:122;1:132;1:133;1:140;1:142;1:158;1:168;1:180;1:201;1:203;1:205;1:208;1:210;1:217;1:228;1:235;1:239;1:245;1:249;1:260;1:266;1:280;1:284;1:288;1:298;1:306;1:310;1:316;1:320;1:327;1:329;1:333;1:335;1:338;1:357;1:368;1:377;	63	0.833333333333333	0.201058201058201	0.251245561590388	0.182321556793955	2.87962538867117	0.00694643445162294	0.0068059	0.967032642054946:4:1;0.650007178987401:3:1;0.105215450692312:2:8;0.00738995402234366:1:53;0.000471151339347942:0:315;
			ENSG00000181935	ENST00000314634	OR4C16	11	6	91	310	6:31;6:135;3:87;3:120;3:194;2:96;2:99;2:112;2:150;2:163;2:186;2:208;2:220;2:264;2:294;1:5;1:9;1:22;1:59;1:61;1:72;1:73;1:90;1:94;1:97;1:98;1:101;1:109;1:116;1:118;1:122;1:126;1:130;1:137;1:138;1:139;1:141;1:143;1:146;1:147;1:149;1:164;1:182;1:190;1:191;1:196;1:214;1:226;1:236;1:241;1:247;1:250;1:253;1:254;1:255;1:257;1:259;1:268;1:274;1:276;1:277;1:284;1:289;1:302;1:304;	65	0.790322580645161	0.293548387096774	0.525200960434282	0.235314086934465	4.27148875371664	0.0144280921838143	0.000000027955	0.999891923117087:6:2;0.607346043154968:3:3;0.0785200850205937:2:10;0.00467228840825762:1:50;0.000258535337030907:0:245;
		Driver_mutation_sites:
			Gene_id	Transcript_id	Chr	Gene_name	Protein_position	z	Q(z)	Protein_mutation
			ENSG00000106278	ENST00000393386	7	PTPRZ1	833	5	0.998756375217605	R833C:4;R833H:1;
			ENSG00000106278	ENST00000393386	7	PTPRZ1	636	4	0.968955612892411	S636L:3;S636T:1;
			ENSG00000106278	ENST00000393386	7	PTPRZ1	205	3	0.548130755870424	R205H:2;R205C:1;
			ENSG00000106278	ENST00000393386	7	PTPRZ1	515	3	0.548130755870424	D515N:3;
			ENSG00000106278	ENST00000393386	7	PTPRZ1	612	3	0.548130755870424	E612K:3;
			ENSG00000105641	ENST00000222248	19	SLC5A5	569	5	0.999903970051843	R569W:5;
			ENSG00000105641	ENST00000222248	19	SLC5A5	438	3	0.896681963701284	P438S:2;P438L:1;		
		The meaning of some parameters in the result:
			Max_hit：the maximum number of mutations at each site of a gene		Count：the total number of mutations in a gene
			Distribution：mutation distribution at all sites of a gene			Mut_sites：all mutation sites in a gene
			f0：the ratio of non-mutated sites in a gene						mean：the average number of mutations at all sites of a gene
			m0：the average number of mutations on the passenger_site			m1：the average number of mutations on the driver_site
			eta：the number of driver_site at all sites of a gene				LRT_p-value：Likelihood Ratio Test value
			Q(z)：posterior probability											z：The number of mutations in the locus	
				