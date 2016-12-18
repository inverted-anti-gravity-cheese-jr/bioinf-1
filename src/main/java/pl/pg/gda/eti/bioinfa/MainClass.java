/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package pl.pg.gda.eti.bioinfa;

import java.awt.BorderLayout;
import javax.swing.JFrame;
import javax.swing.JLabel;

/**
 *
 * @author wojtek
 */
public class MainClass {
    
    public static void main(String[] args) {
	
	JFrame frame = new JFrame();
	
	JLabel label = new JLabel("Hello World!");
	
	frame.setLayout(new BorderLayout());
	frame.add(label, BorderLayout.CENTER);
	frame.setSize(400, 300);
	frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
	frame.setLocationRelativeTo(null);
	frame.setVisible(true);
	
    }
    
}
