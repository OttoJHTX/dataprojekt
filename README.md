# Problemstilling

Molekylærbiologiens centrale dogme beskriver, hvordan det genetiske informationsflow løber fra DNA til RNA til protein. Der findes over 10.000 gener spredt ud i afgrænsede regioner i kromosomerne, som består af DNA. Hvert enkelt gen er forskriften for RNA, som produceres via transskriptionsprocessen. RNA bliver derefter aflæst via translation, hvilket fører til dannelsen af proteiner, som udfører det væld af processer, der foregår i en celle. <br> <br>
<p>
    <img width="450" alt="central_d" src="https://github.com/OttoJHTX/dataprojekt/assets/49984447/69349634-7729-42d3-898b-f45b653eb80e">
    <br>
    <em>Figur 1</em>
</p>
Hvert enkelt trin i det centrale dogme er fundamentalt for alt cellulært liv. Transskriptionen – dannelsen af RNA ved aflæsning af DNA’et fra et givet gen – er i sig selv en kompliceret proces, der består af en lang række koordinerede trin, der udføres af et stort antal proteiner. Transskriptionsprocessen kan deles op i tre faser – initiering, elongering og terminering. I forbindelse med initieringen rekrutteres en polymerase til genet vha. en promoter og påbegynder elongeringen, hvorved genet aflæses i en bestemt retning (se figur 2). Når polymerasen møder en terminator, igangsættes termineringen. Mens initieringen typisk sker fra et defineret område ved promoteren, er termineringen en mere stokastisk proces, hvor polymerasen bremses ned efter mødet med terminatoren og stopper mere eller mindre tilfældigt i nedstrømsregionen. På trods af denne stokasticitet er det dog klart, at processen skal være under stram kontrol, da polymerasen skal være fuldstændigt standset, inden det næste gen begynder. <br> <br>
<p>
    <img width="450" alt="transkription" src="https://github.com/OttoJHTX/dataprojekt/assets/49984447/343be33a-423c-4277-bba3-be89cdda21c7">
    <br>
    <em>Figur 2</em>
</p>
Således er transskriptions-terminering også en kompliceret proces, der involverer over 50 forskellige proteiner. Defekter i termineringen, hvor polymerasen ikke bremses ordentligt, fører til såkaldt “read-through”, hvor transskriptionen fortsætter længere end normalt i den nedstrøms region. Sådanne defekter kan opstå hos patienter, hvor et eller flere af de involverede proteiner er dysfunktionelle. Det kan også opstå i forbindelse med cellulært stress, som ofte opstår midlertidigt i normale individer, men mere ukontrollerbart i forbindelse med sygdomme som for eksempel kræft. Et åbent spørgsmål i feltet er, hvorvidt read-through blot er en konsekvens af stress, eller om det tjener en funktion og derfor er et forsøg fra cellen på at modvirke stress-tilstanden.
For at kunne forstå termineringen og dens involverede trin er det relevant at kunne udføre en detaljeret beskrivelse af read-through under forskellige betingelser. Vi har i dette projekt været med til at designe en model og algoritme, der kan lave en kvantitativ beskrivelse af read-through på individuelle gener baseret på transskriptionsdata fra celler udsat for forskellige betingelser.
<br> <br>

# Data

## Datastruktur
Hele processen for vores data er vist i figur 3. Efter at DNA- eller RNA-prøver er blevet sekventeret i små dele kaldet "reads", bliver de mange millioner sekvenslæs-
ninger opbevaret i en FASTQ-fil. FASTQ-filen er et tekstbaseret filformat, der in-
deholder reads som tekststrenge, eksempelvis sekvensen af nukleotider: "GATTTG
GGGTTC....". Derudover inkluderer hver FASTQ-fil også information om kvaliteten
eller pålideligheden af læsningen for hver enkelt base i sekvensen.
Dette bliver konverteret til en BigWig-fil, som indeholder en værdi for, hvor meget
hver position af genets sekvens er afdækket af reads fra FASTQ-filen, dvs. hvor
meget data der er læst på hver position. Vi har også en annoteringsfil, der giver in-
formation om, hvor i BigWig-filen, de specifikke gener og transskripter er. Vi bruger
så BigWig-filen og annoteringsfilen i sammenhold, så vi kan modellere afdæknin-
gen af specifikke gener og transskripter. <br> <br>
Et eksempel på vores data ses i figur 4, hvor x-aksen er position på genet, og y-
aksen er mængden af data normaliseret. Det ses så hvordan vores kontrol-data, den
røde kurve, har minimal udsving efter terminering (kaldet "downstream of gene"),
mens test-data, den blå kurve, udsvinger meget efter termineringen. Dermed må
der have været en defekt i den blå kurves terminering. <br> <br>
Der er muligvis biases i data, da de maskiner, der læser sekvenserne, kan produc-
ere støj. Derudover er der altid mulighed for menneskelige fejl. Vi vil dog antage, at
biases er minimale, da sekvenserne er blevet læst på et laboratorie i et kontrolleret
miljø.
<p>
    <img width="600" alt="Flow" src="https://github.com/OttoJHTX/dataprojekt/assets/49984447/e982dae5-1eaf-44ae-b24c-c10fd0ef467f">
    <br>
    <em>Figur 3</em>
</p>

## Preprocessering
For at generelisere generne så de passer de er generiske når vi skal modellere på dem, har vi brug for at normalisere vores data. Det første vi gør er $log_2$-transformere data. Her ligger vi først 1 til alle observationer, så når covereage er 0, forbliver det sådan.
<br> <br>
Vi fortsatte derefter til at finde kroppen af generne, så vi kunne undersøge antagelsen om, at kroppen ville være ens for gener med- og uden defekt i termineringen. Udfordringen ved dette var, at vi i annoteringsfilen har flere annotationer af transkripter, hvor der var forskellige koordinater for begyndelsen (TSS) og slutningen (TES) af kroppen. Disse var vi nødt til at sammenligne så vi var sikre på at finde de rigtige. Siden DNA-strenge både kan læses forlæns og baglæns, var vi også nødt til at tage i betragtning hvilken retning vores datapunkter var.



 - Binning
 - Identificér kroppen
 - Antagelser for at kroppen er ens
 - Normaliseret
     - log-transformer, plus 1

# Modellering

## Curve fitting
 - Experiments, gaussian, exponential
 - Plots
 - Double Sigmoid

## Hidden Markov Model
 
# Perspektivering og videreudvikling

 - Neural netværk
     - Diskussion om hvorvidt det er brugbart at bruge det
