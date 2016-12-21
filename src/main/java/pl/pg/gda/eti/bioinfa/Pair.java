/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package pl.pg.gda.eti.bioinfa;

/**
 *
 * @author wojtek
 */
public class Pair<T> {
    private T _1;
    private T _2;
    
    public Pair() {
    }
    
    public Pair(T _1, T _2) {
	this._1 = _1;
	this._2 = _2;
    }

    public T get1() {
	return _1;
    }

    public void set1(T _1) {
	this._1 = _1;
    }

    public T get2() {
	return _2;
    }

    public void set2(T _2) {
	this._2 = _2;
    }
    
}
