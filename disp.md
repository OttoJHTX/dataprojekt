## Disposition:
* Introduktion til projektet: Det centrale dogme, transkriptionelle defekter & hvorfor?
* Data: Hvad arbejder vi med? Hvordan ser det ud? Hvordan arbejder vi med det? - Inkluder normalisering
* Modellering: Double-sigmoidal + Hidden Markov.
* Readthrough Analysis: Eksempel.
* Perspektiv: Machine-learning til detection af readthrough?


## Spørgsmål til Søren:
* Hvad er formålet? Hvad vil du have ud af at have en generel rt-analyse funktion?
  * Grunden er at vi interesserer os for transkriptions-terminering. Forsøge at beskrive dette. Matematisk at beskrive hvad processen gør. Grundvidenskabenlig beskrivelse af essentiel process. Dette skal der være styr på for at celler har det godt. Længde og styrke
  * Når celler er stressede (evt. kræft) har nogen en teori om at genet forsøger at modvirke dette i form af terminerings-defekten.
* Gennemgang af "Double-sigmoidal + Hidden Markov"
  * Independent Gaussian fitting, 2 states - 0 = ingen forskel, 1 = forskel (hidden), Viterbi decoding (robust)
* Hvorfor rt_analyse negativ ikke giver samme output metrics som positiv?
* Hvad var problemet med rt_analyse til at starte med?
* Gennemgang af molekylets historie.
* Forklaring af region og coverage.
* Binner vi stadig data?
  

## To-do:
* Lave summary-stats + plots på data fra det endelige dataset.
* Pimp our GitHub.
* Omskrive og flottificere kode (DEMO).
* Plots produceret med egen kode.
* Få styr på HMM & Double-Sigmoidal.
* Splitte rt analysis op

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
