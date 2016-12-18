/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package pl.pg.gda.eti.bioinfa;

import java.awt.BorderLayout;
import java.io.File;
import javax.swing.JFrame;
import javax.swing.JLabel;
import org.biojava.bio.alignment.AlignmentAlgorithm;
import org.biojava.bio.alignment.AlignmentPair;
import org.biojava.bio.alignment.NeedlemanWunsch;
import org.biojava.bio.alignment.SmithWaterman;
import org.biojava.bio.alignment.SubstitutionMatrix;
import org.biojava.bio.seq.DNATools;
import org.biojava.bio.seq.Sequence;
import org.biojava.bio.symbol.AlphabetManager;
import org.biojava.bio.symbol.FiniteAlphabet;

/**
 *
 * @author wojtek
 */
public class MainClass {

    public static void main(String[] args) {
	try {
	    // The alphabet of the sequences. For this example DNA is choosen.
	    FiniteAlphabet alphabet = (FiniteAlphabet) AlphabetManager.alphabetForName("DNA");
	    // Read the substitution matrix file. 
	    // For this example the matrix NUC.4.4 is good.
	    SubstitutionMatrix matrix = new SubstitutionMatrix(alphabet, new File("mat"));
	    // Define the default costs for sequence manipulation for the global alignment.
	    AlignmentAlgorithm aligner = new NeedlemanWunsch(
		    (short) 0, // match
		    (short) 3, // replace
		    (short) 2, // insert
		    (short) 2, // delete
		    (short) 1, // gapExtend
		    matrix // SubstitutionMatrix
	    );
	    
	    Sequence query = DNATools.createDNASequence("ATAAGC", "query");
	    Sequence target = DNATools.createDNASequence("AAAAACG", "target");
	    // Perform an alignment and save the results.
	    AlignmentPair needleAlignmentPair = aligner.pairwiseAlignment(
		    query, // first sequence
		    target // second one
	    );
	    // Print the alignment to the screen
	    
	    System.out.println("Global alignment with Needleman-Wunsch:\n" + needleAlignmentPair.formatOutput());

	    // Perform a local alginment from the sequences with Smith-Waterman. 
	    // Firstly, define the expenses (penalties) for every single operation.
	    aligner = new SmithWaterman(
		    (short) -1, // match
		    (short) 3, // replace 
		    (short) 2, // insert
		    (short) 2, // delete
		    (short) 1, // gapExtend
		    matrix // SubstitutionMatrix
	    );
	    // Perform the local alignment.
	    AlignmentPair smithAlignmentPair = aligner.pairwiseAlignment(query, target);
	    System.out.println("\nlocal alignment with SmithWaterman:\n" + smithAlignmentPair.formatOutput());
	} catch (Exception exc) {
	    exc.printStackTrace();
	}

	JFrame frame = new JFrame();

	JLabel label = new JLabel("Hello World!");

	frame.setLayout(new BorderLayout());
	frame.add(label, BorderLayout.CENTER);
	frame.setSize(400, 300);
	frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
	frame.setLocationRelativeTo(null);
	frame.setVisible(true);

	FiniteAlphabet alphabet = (FiniteAlphabet) AlphabetManager.alphabetForName("DNA");
	System.out.println(alphabet);

    }

}
