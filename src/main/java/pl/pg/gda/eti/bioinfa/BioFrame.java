/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package pl.pg.gda.eti.bioinfa;

import javax.swing.JFileChooser;
import org.biojava.bio.seq.DNATools;
import org.biojava.bio.symbol.SymbolList;

/**
 *
 * @author wojtek
 */
public class BioFrame extends javax.swing.JFrame {

    /**
     * Creates new form BioFrame
     */
    public BioFrame() {
	initComponents();
    }

    /**
     * This method is called from within the constructor to initialize the form.
     * WARNING: Do NOT modify this code. The content of this method is always
     * regenerated by the Form Editor.
     */
    @SuppressWarnings("unchecked")
    // <editor-fold defaultstate="collapsed" desc="Generated Code">//GEN-BEGIN:initComponents
    private void initComponents() {

        jFileChooser1 = new javax.swing.JFileChooser();
        firstSequenceTextLabel = new java.awt.Label();
        firstSequenceTextField = new javax.swing.JTextField();
        secondSequenceTextLabel = new java.awt.Label();
        secondSequenceTextField = new javax.swing.JTextField();
        matrixTextLabel = new java.awt.Label();
        matrixTextField = new javax.swing.JTextField();
        matrixFileChooser = new javax.swing.JButton();
        globalAlignmentButton = new javax.swing.JButton();
        localAlignmentButton = new javax.swing.JButton();
        editDistanceButton = new javax.swing.JButton();
        resultTextLabel = new java.awt.Label();
        jScrollPane1 = new javax.swing.JScrollPane();
        resultTextArea = new javax.swing.JTextArea();
        rnaCheckbox = new javax.swing.JCheckBox();

        setDefaultCloseOperation(javax.swing.WindowConstants.EXIT_ON_CLOSE);

        firstSequenceTextLabel.setText("Sekwencja pierwsza");

        firstSequenceTextField.setText("ATAAGC");
        firstSequenceTextField.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                firstSequenceTextFieldActionPerformed(evt);
            }
        });

        secondSequenceTextLabel.setText("Sekwencja druga");

        secondSequenceTextField.setText("AAAAGC");

        matrixTextLabel.setText("Macierz podobieństwa");

        matrixTextField.setEditable(false);
        matrixTextField.setText("nucl");

        matrixFileChooser.setText("Wybierz plik");
        matrixFileChooser.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                matrixFileChooserActionPerformed(evt);
            }
        });

        globalAlignmentButton.setText("Pokaż dopasowanie globalne");
        globalAlignmentButton.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                globalAlignmentButtonActionPerformed(evt);
            }
        });

        localAlignmentButton.setText("Pokaż dopasowanie lokalne");
        localAlignmentButton.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                localAlignmentButtonActionPerformed(evt);
            }
        });

        editDistanceButton.setText("Pokaż odległość edycyjną");
        editDistanceButton.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                editDistanceButtonActionPerformed(evt);
            }
        });

        resultTextLabel.setText("Wynik");

        resultTextArea.setColumns(20);
        resultTextArea.setRows(5);
        jScrollPane1.setViewportView(resultTextArea);

        rnaCheckbox.setText("Czy jest sekwencją RNA?");
        rnaCheckbox.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                rnaCheckboxActionPerformed(evt);
            }
        });

        javax.swing.GroupLayout layout = new javax.swing.GroupLayout(getContentPane());
        getContentPane().setLayout(layout);
        layout.setHorizontalGroup(
            layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(layout.createSequentialGroup()
                .addContainerGap()
                .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    .addComponent(jScrollPane1, javax.swing.GroupLayout.Alignment.TRAILING)
                    .addComponent(globalAlignmentButton, javax.swing.GroupLayout.DEFAULT_SIZE, 473, Short.MAX_VALUE)
                    .addComponent(firstSequenceTextField)
                    .addComponent(secondSequenceTextField)
                    .addGroup(layout.createSequentialGroup()
                        .addComponent(matrixTextField)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addComponent(matrixFileChooser))
                    .addComponent(localAlignmentButton, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                    .addComponent(editDistanceButton, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                    .addGroup(layout.createSequentialGroup()
                        .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                            .addComponent(firstSequenceTextLabel, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                            .addComponent(secondSequenceTextLabel, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                            .addComponent(matrixTextLabel, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                            .addComponent(resultTextLabel, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))
                        .addGap(0, 0, Short.MAX_VALUE))
                    .addComponent(rnaCheckbox, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE))
                .addContainerGap())
        );
        layout.setVerticalGroup(
            layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(layout.createSequentialGroup()
                .addGap(15, 15, 15)
                .addComponent(firstSequenceTextLabel, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(firstSequenceTextField, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(secondSequenceTextLabel, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(secondSequenceTextField, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(matrixTextLabel, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(matrixTextField, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(matrixFileChooser))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(globalAlignmentButton)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(localAlignmentButton)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(editDistanceButton)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(rnaCheckbox)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(resultTextLabel, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(jScrollPane1)
                .addContainerGap())
        );

        pack();
    }// </editor-fold>//GEN-END:initComponents

    private void globalAlignmentButtonActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_globalAlignmentButtonActionPerformed

	try {
	    String result = "";
	    String sequence1 = firstSequenceTextField.getText();
	    String sequence2 = secondSequenceTextField.getText();
	    String similarityMatrixFile = matrixTextField.getText();
	    Pair<SymbolList> globalAlignment;
	    if (rnaCheckbox.isSelected()) {
		Pair<Pair<SymbolList>> alignmentWithProtein = BioAlgorithms.getGlobalAlignmentForRNA(sequence1, sequence2, similarityMatrixFile);
		result += "Aminokwasy:\n" + alignmentWithProtein.get2().get1().seqString() + "\n" + alignmentWithProtein.get2().get2().seqString() + "\n\n";
		globalAlignment = alignmentWithProtein.get1();
	    } else {
		globalAlignment = BioAlgorithms.getGlobalAlignment(sequence1, sequence2, similarityMatrixFile);
	    }
	    result += "Wynik:\n" + globalAlignment.get1().seqString() + "\n" + globalAlignment.get2().seqString();
	    resultTextArea.setText(result);
	} catch (Exception e) {
	    e.printStackTrace();
	}
    }//GEN-LAST:event_globalAlignmentButtonActionPerformed

    private void localAlignmentButtonActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_localAlignmentButtonActionPerformed
	try {
	    String result = "";
	    String sequence1 = firstSequenceTextField.getText();
	    String sequence2 = secondSequenceTextField.getText();
	    String similarityMatrixFile = matrixTextField.getText();
	    Pair<SymbolList> localAlignment;
	    if (rnaCheckbox.isSelected()) {
		Pair<Pair<SymbolList>> alignmentWithProtein = BioAlgorithms.getLocalAlignmentForRNA(sequence1, sequence2, similarityMatrixFile);
		result += "Aminokwasy:\n" + alignmentWithProtein.get2().get1().seqString() + "\n" + alignmentWithProtein.get2().get2().seqString() + "\n\n";
		localAlignment = alignmentWithProtein.get1();
	    }
	    else {
		localAlignment = BioAlgorithms.getLocalAlignment(sequence1, sequence2, similarityMatrixFile);
	    }
	    result += "Wynik:\n" + localAlignment.get1().seqString() + "\n" + localAlignment.get2().seqString();
	    resultTextArea.setText(result);
	} catch (Exception e) {
	    e.printStackTrace();
	}
    }//GEN-LAST:event_localAlignmentButtonActionPerformed

    private void matrixFileChooserActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_matrixFileChooserActionPerformed
	try {
	    JFileChooser chooser = new JFileChooser();
	    chooser.showOpenDialog(this);
	    matrixTextField.setText(chooser.getSelectedFile().getAbsolutePath());
	}
	catch (Exception e) {
	}
    }//GEN-LAST:event_matrixFileChooserActionPerformed

    private void firstSequenceTextFieldActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_firstSequenceTextFieldActionPerformed
	// TODO add your handling code here:
    }//GEN-LAST:event_firstSequenceTextFieldActionPerformed

    private String lastDnaSeq1;
    private String lastDnaSeq2;

    private void rnaCheckboxActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_rnaCheckboxActionPerformed
	if (rnaCheckbox.isSelected()) {
	    try {
		lastDnaSeq1 = firstSequenceTextField.getText();
		lastDnaSeq2 = secondSequenceTextField.getText();
		SymbolList rna = DNATools.toRNA(DNATools.createDNA(firstSequenceTextField.getText()));
		firstSequenceTextField.setText(rna.seqString().toUpperCase());
		rna = DNATools.toRNA(DNATools.createDNA(secondSequenceTextField.getText()));
		secondSequenceTextField.setText(rna.seqString().toUpperCase());
	    } catch (Exception ex) {
		ex.printStackTrace();
	    }
	} else {
	    firstSequenceTextField.setText(lastDnaSeq1);
	    secondSequenceTextField.setText(lastDnaSeq2);
	}
    }//GEN-LAST:event_rnaCheckboxActionPerformed

    private void editDistanceButtonActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_editDistanceButtonActionPerformed
        try {
	    String sequence1 = firstSequenceTextField.getText();
	    String sequence2 = secondSequenceTextField.getText();
	    String similarityMatrixFile = matrixTextField.getText();
	    if(rnaCheckbox.isSelected()) {
		int editDistance = BioAlgorithms.getEditDistanceForRNA(sequence1, sequence2, similarityMatrixFile);
		resultTextArea.setText("Odległość edycyjna:\n" + editDistance);
	    }
	    else {
		int editDistance = BioAlgorithms.getEditDistance(sequence1, sequence2, similarityMatrixFile);
		resultTextArea.setText("Odległość edycyjna:\n" + editDistance);
	    }
	} catch (Exception e) {
	}
    }//GEN-LAST:event_editDistanceButtonActionPerformed


    // Variables declaration - do not modify//GEN-BEGIN:variables
    private javax.swing.JButton editDistanceButton;
    private javax.swing.JTextField firstSequenceTextField;
    private java.awt.Label firstSequenceTextLabel;
    private javax.swing.JButton globalAlignmentButton;
    private javax.swing.JFileChooser jFileChooser1;
    private javax.swing.JScrollPane jScrollPane1;
    private javax.swing.JButton localAlignmentButton;
    private javax.swing.JButton matrixFileChooser;
    private javax.swing.JTextField matrixTextField;
    private java.awt.Label matrixTextLabel;
    private javax.swing.JTextArea resultTextArea;
    private java.awt.Label resultTextLabel;
    private javax.swing.JCheckBox rnaCheckbox;
    private javax.swing.JTextField secondSequenceTextField;
    private java.awt.Label secondSequenceTextLabel;
    // End of variables declaration//GEN-END:variables
}
