/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package pl.pg.gda.eti.bioinfa;

import java.io.File;
import java.util.Set;
import org.biojava.bio.alignment.AlignmentAlgorithm;
import org.biojava.bio.alignment.AlignmentPair;
import org.biojava.bio.alignment.NeedlemanWunsch;
import org.biojava.bio.alignment.SmithWaterman;
import org.biojava.bio.alignment.SubstitutionMatrix;
import org.biojava.bio.seq.DNATools;
import org.biojava.bio.seq.RNATools;
import org.biojava.bio.seq.Sequence;
import org.biojava.bio.symbol.AlphabetManager;
import org.biojava.bio.symbol.FiniteAlphabet;
import org.biojava.bio.symbol.IllegalSymbolException;
import org.biojava.bio.symbol.ManyToOneTranslationTable;
import org.biojava.bio.symbol.SimpleAtomicSymbol;
import org.biojava.bio.symbol.Symbol;
import org.biojava.bio.symbol.SymbolList;

/**
 *
 * @author wojtek
 */
public class BioAlgorithms {

    public static Pair<SymbolList> getLocalAlignment(String sequence1, String sequence2, String similarityMatrixFile) throws Exception {
	// The alphabet of the sequences. For this example DNA is choosen.
	FiniteAlphabet alphabet = (FiniteAlphabet) AlphabetManager.alphabetForName("DNA");
	// Read the substitution matrix file. 
	SubstitutionMatrix matrix = new SubstitutionMatrix(alphabet, new File(similarityMatrixFile));
	// Define the default costs for sequence manipulation for the global alignment.
	short noCost = 0;
	AlignmentAlgorithm aligner = new SmithWaterman(noCost, noCost, noCost, noCost, noCost, matrix);
	
	Sequence query = DNATools.createDNASequence(sequence1, "query");
	Sequence target = DNATools.createDNASequence(sequence2, "target");
	
	// Perform an alignment and save the results.
	AlignmentPair needleAlignmentPair = aligner.pairwiseAlignment(query, target);
	
	return new Pair<SymbolList>(needleAlignmentPair.getQuery(), needleAlignmentPair.getSubject());
    }
    
    public static Pair<SymbolList> getGlobalAlignment(String sequence1, String sequence2, String similarityMatrixFile) throws Exception {
	// The alphabet of the sequences. For this example DNA is choosen.
	FiniteAlphabet alphabet = (FiniteAlphabet) AlphabetManager.alphabetForName("DNA");
	// Read the substitution matrix file. 
	SubstitutionMatrix matrix = new SubstitutionMatrix(alphabet, new File(similarityMatrixFile));
	// Define the default costs for sequence manipulation for the global alignment.
	short noCost = 0;
	AlignmentAlgorithm aligner = new NeedlemanWunsch(noCost, noCost, noCost, noCost, noCost, matrix);

	Sequence query = DNATools.createDNASequence(sequence1, "query");
	Sequence target = DNATools.createDNASequence(sequence2, "target");
	// Perform an alignment and save the results.
	AlignmentPair needleAlignmentPair = aligner.pairwiseAlignment(query, target);
	
	return new Pair<SymbolList>(needleAlignmentPair.getQuery(), needleAlignmentPair.getSubject());
    }
    
    public static Pair<Pair<SymbolList>> getGlobalAlignmentForRNA(String sequence1, String sequence2, String similarityMatrixFile) throws Exception {
	// The alphabet of the sequences. For this example DNA is choosen.
	FiniteAlphabet alphabet = (FiniteAlphabet) AlphabetManager.alphabetForName("PROTEIN-TERM");
	// Read the substitution matrix file. 
	SubstitutionMatrix matrix = new SubstitutionMatrix(alphabet, new File(similarityMatrixFile));
	// Define the default costs for sequence manipulation for the global alignment.
	short noCost = 0;
	AlignmentAlgorithm aligner = new NeedlemanWunsch(noCost, noCost, noCost, noCost, noCost, matrix);

	Sequence query = RNATools.createRNASequence(sequence1, "query");
	Sequence target = RNATools.createRNASequence(sequence2, "target");
	
	SymbolList queryProtein = RNATools.translate(query);
	SymbolList targetProtein = RNATools.translate(target);
	
	// Perform an alignment and save the results.
	AlignmentPair needleAlignmentPair = aligner.pairwiseAlignment(queryProtein, targetProtein);
	
	SymbolList alignedQuery = untranslateProtein(query, needleAlignmentPair.getQuery());
	SymbolList alignedSubject = untranslateProtein(query, needleAlignmentPair.getSubject());
	
	Pair<SymbolList> pairResult = new Pair<SymbolList>(alignedQuery, alignedSubject);
	Pair<SymbolList> pairProtein = new Pair<SymbolList>(queryProtein, targetProtein);
	
	return new Pair<Pair<SymbolList>>(pairResult, pairProtein);
    }
    
    public static Pair<Pair<SymbolList>> getLocalAlignmentForRNA(String sequence1, String sequence2, String similarityMatrixFile) throws Exception {
	// The alphabet of the sequences. For this example DNA is choosen.
	FiniteAlphabet alphabet = (FiniteAlphabet) AlphabetManager.alphabetForName("PROTEIN-TERM");
	// Read the substitution matrix file. 
	SubstitutionMatrix matrix = new SubstitutionMatrix(alphabet, new File(similarityMatrixFile));
	// Define the default costs for sequence manipulation for the global alignment.
	short noCost = 0;
	AlignmentAlgorithm aligner = new SmithWaterman(noCost, noCost, noCost, noCost, noCost, matrix);

	Sequence query = RNATools.createRNASequence(sequence1, "query");
	Sequence target = RNATools.createRNASequence(sequence2, "target");
	
	SymbolList queryProtein = RNATools.translate(query);
	SymbolList targetProtein = RNATools.translate(target);
	
	// Perform an alignment and save the results.
	AlignmentPair needleAlignmentPair = aligner.pairwiseAlignment(queryProtein, targetProtein);
	
	SymbolList alignedQuery = untranslateProtein(query, needleAlignmentPair.getQuery());
	SymbolList alignedSubject = untranslateProtein(query, needleAlignmentPair.getSubject());
	
	Pair<SymbolList> pairResult = new Pair<SymbolList>(alignedQuery, alignedSubject);
	Pair<SymbolList> pairProtein = new Pair<SymbolList>(queryProtein, targetProtein);
	
	return new Pair<Pair<SymbolList>>(pairResult, pairProtein);
    }
       
    public static int getEditDistance(String sequence1, String sequence2, String similarityMatrixFile) throws Exception {
        // The alphabet of the sequences. For this example DNA is choosen.
	FiniteAlphabet alphabet = (FiniteAlphabet) AlphabetManager.alphabetForName("DNA");
	// Read the substitution matrix file. 
	SubstitutionMatrix matrix = new SubstitutionMatrix(alphabet, new File(similarityMatrixFile));
	// Define the default costs for sequence manipulation for the global alignment.
	short noCost = 0;
	NeedlemanWunsch aligner = new NeedlemanWunsch(noCost, noCost, noCost, noCost, noCost, matrix);

	Sequence query = DNATools.createDNASequence(sequence1, "query");
	Sequence target = DNATools.createDNASequence(sequence2, "target");
	// Perform an alignment and save the results.
	AlignmentPair needleAlignmentPair = aligner.pairwiseAlignment(query, target);
        
        int editDistance = needleAlignmentPair.length() - needleAlignmentPair.getNumIdenticals();
	return editDistance;
    }
    
    public static int getEditDistanceForRNA(String sequence1, String sequence2, String similarityMatrixFile) throws Exception {
        // The alphabet of the sequences. For this example DNA is choosen.
	FiniteAlphabet alphabet = (FiniteAlphabet) AlphabetManager.alphabetForName("PROTEIN-TERM");
	// Read the substitution matrix file. 
	SubstitutionMatrix matrix = new SubstitutionMatrix(alphabet, new File(similarityMatrixFile));
	// Define the default costs for sequence manipulation for the global alignment.
	short noCost = 0;
	AlignmentAlgorithm aligner = new NeedlemanWunsch(noCost, noCost, noCost, noCost, noCost, matrix);

	Sequence query = RNATools.createRNASequence(sequence1, "query");
	Sequence target = RNATools.createRNASequence(sequence2, "target");
	
	SymbolList queryProtein = RNATools.translate(query);
	SymbolList targetProtein = RNATools.translate(target);
	
	// Perform an alignment and save the results.
	AlignmentPair needleAlignmentPair = aligner.pairwiseAlignment(queryProtein, targetProtein);
        
        int editDistance = needleAlignmentPair.length() - needleAlignmentPair.getNumIdenticals();
	return editDistance;
    }
    
    public static SymbolList untranslateProtein(SymbolList oldRna, SymbolList protein) throws IllegalSymbolException {
	
	// inicjalizacja
	ManyToOneTranslationTable transcriptionTable = RNATools.getGeneticCode("UNIVERSAL");
	String rnaHelper = "";
	
	int searchFrom = 0;
	int searchTo = 0;
	
	try {
	    // dla każdego znaku w aminkwasie...
	    for(int i = 1; i <= protein.length(); i++) {
		// przetłumacz znak na RNA (będzie kilka kombinacji)
		Set symbolSet = transcriptionTable.untranslate(protein.symbolAt(i));
		
		// jeśli pusty to znaczy że usunięty lub dodany (czyli -)
		if(symbolSet.isEmpty()) {
		    rnaHelper += "-";
		    searchTo ++;
		}
		// jeśli są jakieś kombinacje to wybierz właściwą
		else {
		    String correctSymbol = null;
		    int minInd = searchTo + 1;
		    // dla wszystkich kombinacji...
		    for(Object symbolSubsetObject : symbolSet) {
			// konwersja na string
			String symbolStr = symbolToString((SimpleAtomicSymbol) symbolSubsetObject);
			// jak nie znajdzie to weź dowolny symbol
			if(correctSymbol == null) {
			    correctSymbol = symbolStr;
			}
			
			// przeszukaj fragment sekwencji RNA pod kątem występowania symbolu
			int ind = searchInRnaSubstring(symbolStr, oldRna, searchFrom, searchTo);
			
			// wybierz ten który był znaleziony "najwcześniej" w sekwencji
			if (ind < minInd && ind >= 0) {
			    minInd = ind;
			    correctSymbol = symbolStr;
			}
		    }
		    
		    // wstaw nowy symbol (dowolny lub konkretny)
		    rnaHelper += correctSymbol;
		    // jeżeli znaleziono w sekwencji to przesuń wyszukiwanie dla następnych
		    if(minInd < searchTo + 1) {
			System.out.println("Znaleziono symbol '" + correctSymbol + "' na indeksie " + minInd);
			searchFrom = minInd;
			searchTo = Math.max(searchTo, minInd);
		    }
		}
	    }
	    SymbolList newList = RNATools.createRNA(rnaHelper);
	    return newList;
	}
	catch (Exception e) {
	    e.printStackTrace();
	}
	return null;
    }
    
    private static int searchInRnaSubstring(String symbol, SymbolList rna, int searchFrom, int searchTo) {
	// przeszukaj fragment sekwencji RNA pod kątem występowania symbolu
	String rnaString = rna.seqString();
	for(int i = searchFrom; i <= searchTo; i++) {
	    String rnaSub = rnaString.substring(i * 3, (i + 1) * 3);
	    if(symbol.toUpperCase().equals(rnaSub.toUpperCase())) {
		return i;
	    }
	}
	return -1;
    }
    
    private static String symbolToString(SimpleAtomicSymbol sym) {
	// dla symbolu pobierz jego wsyzstkie "podsymbole" i pobierz pierwsze litery pełnich nazw
	// w UPPERCASE
	String rnaSymbol = "";
	for(Object symbolObj: sym.getSymbols()) {
	    Symbol symbol = (Symbol) symbolObj;
	    rnaSymbol += symbol.getName().toUpperCase().charAt(0);
	}
	return rnaSymbol;
    }
}