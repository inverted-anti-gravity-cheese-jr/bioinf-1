/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package pl.pg.gda.eti.bioinfa;

import java.io.File;
import org.biojava.bio.BioException;
import org.biojava.bio.alignment.AlignmentAlgorithm;
import org.biojava.bio.alignment.AlignmentPair;
import org.biojava.bio.alignment.NeedlemanWunsch;
import org.biojava.bio.alignment.SmithWaterman;
import org.biojava.bio.alignment.SubstitutionMatrix;
import org.biojava.bio.seq.DNATools;
import org.biojava.bio.seq.ProteinTools;
import org.biojava.bio.seq.RNATools;
import org.biojava.bio.seq.Sequence;
import org.biojava.bio.symbol.AlphabetManager;
import org.biojava.bio.symbol.FiniteAlphabet;
import org.biojava.bio.symbol.IllegalSymbolException;
import org.biojava.bio.symbol.SimpleReversibleTranslationTable;
import org.biojava.bio.symbol.SymbolList;
import org.biojava.bio.symbol.SymbolListViews;

/**
 *
 * @author wojtek
 */
public class BioAlgorithms {

    public static AlignmentPair getLocalAlignment(String sequence1, String sequence2, String similarityMatrixFile) throws Exception {
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
	
	System.out.println(needleAlignmentPair.formatOutput());
	
	return needleAlignmentPair;
    }
    
    public static AlignmentPair getGlobalAlignment(String sequence1, String sequence2, String similarityMatrixFile) throws Exception {
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
	
	System.out.println(needleAlignmentPair.formatOutput());

	return needleAlignmentPair;
    }
    
    public static AlignmentPair getGlobalAlignmentForRNA(String sequence1, String sequence2, String similarityMatrixFile) throws Exception {
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
	
	System.out.println(needleAlignmentPair.formatOutput());
	
	return needleAlignmentPair;
    }
    
    public static AlignmentPair getLocalAlignmentForRNA(String sequence1, String sequence2, String similarityMatrixFile) throws Exception {
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
	
	System.out.println(needleAlignmentPair.formatOutput());
	
	return needleAlignmentPair;
    }

}
