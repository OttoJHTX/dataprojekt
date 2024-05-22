# Problemstilling

Det centrale dogme inden for biologi siger, at genetisk information går fra DNA
til RNA gennem transskription, og så til protein gennem translation, hvilket kan
ses i figur 1. Transskriptions processen (afbildet i figur 2.) består af initiering,
elongering og terminering, som henholdsvis markerer indledningen af transskrip-
tion, opbygning af RNA-molekylet og afslutningen af transskriptionen. I visse til-
fælde, blandt andet forårsaget af stress, kan der ske defekter i termineringen, så
transskriptionen fortsætter for længe. Vores projekt går ud på at lave en model, der
kan identificere og beskrive disse defekter.

<p>
    <img width="450" alt="central_d" src="https://github.com/OttoJHTX/dataprojekt/assets/49984447/69349634-7729-42d3-898b-f45b653eb80e">
    <br>
    <em>Figur 1</em>
</p>


<p>
    <img width="450" alt="transkription" src="https://github.com/OttoJHTX/dataprojekt/assets/49984447/343be33a-423c-4277-bba3-be89cdda21c7">
    <br>
    <em>Figur 2</em>
</p>

Der eksisterer få værktøjer til analyse af terminerings-defekter i transskription, og
de få, der eksisterer, er mangelfulde og begrænsede i deres evner til at kvantificere
og beskrive defekten. For at adressere denne mangel, vil vores projekt fokusere på
udviklingen af et værktøj der kan modellere og beskrive transskriptionsterminerin-
gen, med henblik på at afgøre om der er en defekt eller ej, såvel som at kvantificere
og beskrive denne.

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
    <em>Figur 2</em>
</p>

## Preprocessering
For at generelisere generne så de passer de er generiske når vi skal modellere på dem, har vi brug for at normalisere vores data. Det første vi gør er $log_2$-transformere data. Her ligger vi først 1 til alle observationer, så når covereage er 0, forbliver det sådan.
<br> <br>
Vi fortsatte derefter til at finde kroppen af generne, så vi kunne undersøge antagelsen om, at kroppen ville være ens for gener med- og uden defekt i termineringen. Udfordringen ved dette var, at vi i annoteringsfilen har flere annoteringer af transkripter, hvor der var forskellige koordinater for begyndelsen (TSS) og slutningen (TES) af kroppen. Disse var vi nødt til at sammenligne så vi var sikre på at finde de rigtige. 
```{pseudo}
if (estimate_TES and transcript_annotations.size() > 1) {
    Determine outer boundaries for TES and TSS estimation
    Extract unique start and end coordinates from transcript annotations
    Initialize vectors for storing differences in expression levels
    Iterate over transcript start coordinates:
        Calculate expression differences for nearby regions
        Store minimum difference in a vector
    Iterate over transcript end coordinates:
        Calculate expression differences for nearby regions
        Store minimum difference in a vector
    Find indices of maximum difference for TES and TSS determination
    
    Update final transcript range:
        Extract TSS and TES coordinates based on maximum difference indices
        Adjust final range accordingly
    
    Determine upstream and downstream TSS and TES coordinates based on strand
    Retrieve corresponding data rows for TSS, TES, upstream TSS, and downstream TES
}
```



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
