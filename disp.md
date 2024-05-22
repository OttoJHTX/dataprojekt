## Disposition:
* Introduktion til projektet: Det centrale dogme, transkriptionelle defekter & hvorfor?
* Data: Hvad arbejder vi med? Hvordan ser det ud? Hvordan arbejder vi med det? - Inkluder normalisering
* Modellering: Double-sigmoidal + Hidden Markov.
* Readthrough Analysis: Eksempel.
* Perspektiv: Machine-learning til detection af readthrough?


## Spørgsmål til Søren:
* Hvad er formålet? Hvad vil du have ud af at have en generel rt-analyse funktion?
* Gennemgang af "Double-sigmoidal + Hidden Markov"
* Hvorfor rt_analyse negativ ikke giver samme output metrics som positiv?
* Hvad var problemet med rt_analyse til at starte med?


## Hvad har vi lavet? 
1. Læste op på det centrale dogme/genetisk struktur generelt.
2. Forståelse af data (metrics), filtyper etc.
3. Forbinde til cluster/server med SSH og lære at navigere derinde/benytte funktionaliteter.
4. Løsning af små opgaver, vi fik fra Søren med brug af data (både coverage og annotering) i R vha. div. BioConductor pakker.
5. Sideløbende lavet vores projektbeskrivelse.
6. Fitte modeller til data (double-sigmoidal, eksponentiel, gaussian) og evaluere.
7. Vi forsøgte at binde de forskellige koncepter sammen (scripts/kode etc. vi har arbejdet på indtil videre) til en funktion der var 'plug & play'.
8. Diskussion med Søren omkring værdien ved at fortsætte med at udarbejde et RT-analyse script, når hans eget RT-analyse script i forvejen er så 'udviklet'. Her bliver ideen om, at introducere en machine-learning dimension introduceret.
9. Vi aftaler at tage projektet i en anden retning - Vi vil benytte output (statistikker, vigtige mål etc.) fra Sørens RT-analyse som input til et neuralt netværk, der skulle have til formål at afgøre om der er RT eller ej.
10. Vi forbereder kode til at lave, teste og gemme netværk vha. PyTorch & Keras.
11. I samarbejde med Søren, videreudvikler vi hans rt_analysis script - vi udvikler både tekniske aspekter samt generaliserer scriptes 'kørsel', så det er mere fleksibelt og kan håndtere at køre n-gener på samme tid.
12. Vi laver også en wrapper til rt_analyse, så det kan køres på clusteret.
13. Vi kører rt_analyse på ca. 10.386 gener (både en med negativ og en med positiv rt for hvert gen) -> får et datasæt ud. 
14. Vi opdager at funktionen, der detekter readthrough ikke måler de nødvendige inputs til netværket hvis der ikke er et readthrough
15. Vi går videre uden et neuralt netværk
