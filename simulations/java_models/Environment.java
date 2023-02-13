// this Java file contains the methods for the actual Environment, and here it initializes the lattice structure of 9 non stem cancer cells and 1 stem cell, and simulated for 60 timesteps

package simulations.java_models;

import java.util.ArrayList;
import java.util.Random;
 
public class Environment {
   public static void main(String[] args) {
       ArrayList<Cell> allCells = new ArrayList<Cell>(); // list of all current cells
       Cell[][] lattice = new Cell[10][10]; // the grid-like lattice of cells
       Random rand = new Random();
      
       // initialize the lattice with 8 non-stem cells and 1 stem-cell
       for (int i = 4; i < 7; i++) {
           for (int j = 4; j < 7; j++) {
               if (i != 5 || j != 5) {
                   lattice[i][j] = new Cell(i, j, false, 10, rand.nextDouble(), rand.nextDouble());
                   allCells.add(lattice[i][j]);
               } else {
                   lattice[i][j] = new Cell(i, j, true, 10, rand.nextDouble(), rand.nextDouble());
                   allCells.add(lattice[i][j]);
               }
           }
       }
      
       // simulate the procedure for 60 time units
       for (int time = 1; time < 60; time++) {
           System.out.println("time = " + time);
           simulate(lattice, allCells); // go through one simulation
           print(lattice); // print out the lattice
           allCells = getAllCells(lattice); // get the new list of current cells in the lattice
       }
   }
  
   /**
    * Returns an ArrayList of all Cells in the current lattice in O(N^2) time.
    * @param lattice: a Cell[][] object that is the current state of the environment
    * @return ArrayList<Cell> of all cells in the lattice
    */
   static ArrayList<Cell> getAllCells(Cell[][] lattice) {
       ArrayList<Cell> allCells = new ArrayList<Cell>();
      
       for (int r = 0; r < lattice.length; r++) {
           for (int c = 0; c < lattice[0].length; c++) {
               if (lattice[r][c] != null) allCells.add(lattice[r][c]);
           }
       }
       return allCells;
   }
  
   /**
    * Simulates through one time-stamp for all current cells in the lattice.
    * @param lattice: a Cell[][] object that is the current state of the environment
    * @param allCells: a ArrayList<Cell> object that has the current list of Cells in the lattice
    */
   static void simulate(Cell[][] lattice, ArrayList<Cell> allCells) {
       // if the list is empty, don't simulate
       if (allCells.isEmpty()) {
           return;
       }
      
       // until the list is empty, do this for every Cell
       while (!allCells.isEmpty()) {
           // choose a random Cell in the list
           Random rand = new Random();
           int length = allCells.size();
           int cellIndex = rand.nextInt(length);
           Cell c = allCells.get(cellIndex);
          
           // probability that the cell randomly dies
           double death = rand.nextDouble();
           if (death < c.spontaneousDeath) {
               // remove the cell from both the list and lattice
               allCells.remove(cellIndex);
               lattice[c.x][c.y] = null;
               continue; // go on to the next cell
           }
          
           // find all free spots around the Cell in the lattice
           ArrayList<ArrayList<Integer>> freeSpots = anyFreeSpot(lattice, c);
           // only if there are free spots, divide the cell
           if (freeSpots.size() > 0) {
               // probability that the cell randomly dies during cell division
               double n = rand.nextDouble();
               if (n < c.spontaneousDeath) {
                   allCells.remove(cellIndex);
                   lattice[c.x][c.y] = null;
                   continue;
               }
              
               // if the cell is NOT a stem cell
               if (!c.isStem) {
                   c.telomereLength--; // reduce its telomere length. if it is less than 0, then it dies.
                   if (c.telomereLength <= 0) {
                       allCells.remove(cellIndex);
                       lattice[c.x][c.y] = null;
                       continue;
                   }
                  
                   // choose a random free spot to put the daughter cell
                   System.out.println("made non-stem daughter cell");
                   int randIndex = rand.nextInt(freeSpots.size());
                   Cell daughterCell = new Cell(freeSpots.get(randIndex).get(0), freeSpots.get(randIndex).get(1),
                           false, c.telomereLength, c.symmetricDivision, c.spontaneousDeath);
                  
                   lattice[freeSpots.get(randIndex).get(0)][freeSpots.get(randIndex).get(0)] = daughterCell;
                   allCells.add(daughterCell); // add daughter cell to the list
               }
               // if the cell IS a stem cell
               else {
                   System.out.println("made stem daughter cell");
                   // get probability if there's a symmetric division
                   double symmetric = rand.nextDouble();
                   int randIndex = rand.nextInt(freeSpots.size());
                   Cell daughterCell;
                   // if the division is symmetric
                   if (symmetric < c.symmetricDivision) { // make two stem daughter cells that are the same
                       daughterCell = new Cell(freeSpots.get(randIndex).get(0), freeSpots.get(randIndex).get(1),
                               true, c.telomereLength, c.symmetricDivision, c.spontaneousDeath);
                   } else { // otherwise make a non-stem daughter cell
                       daughterCell = new Cell(freeSpots.get(randIndex).get(0), freeSpots.get(randIndex).get(1),
                               false, c.telomereLength, c.symmetricDivision, c.spontaneousDeath);
                   }
                  
                   // add daughter cell to free spot
                   lattice[freeSpots.get(randIndex).get(0)][freeSpots.get(randIndex).get(0)] = daughterCell;
                   allCells.add(daughterCell);
               }
           }
          
           // remove it from the list of current cells and go to the next cell
           allCells.remove(cellIndex);
       }
   }
  
   /**
    * Finds and returns the number of free spots around a Cell c in the lattice.
    * @param lattice: a Cell[][] object that is the current state of the environment
    * @param c: the current Cell
    * @return ArrayList<ArrayList<Integer>>: all the positions around a cell that are currently empty.
    */
   static ArrayList<ArrayList<Integer>> anyFreeSpot(Cell[][] lattice, Cell c) {
       // initialize all the returned array
       ArrayList<ArrayList<Integer>> freeSpots = new ArrayList<ArrayList<Integer>>();
       int x = c.x;
       int y = c.y;
      
       // for the eight spots around the cell
       for (int xx = -1; xx <= 1; xx++) {
           for (int yy = -1; yy <= 1; yy++) {
               if (xx == 0 && yy == 0) {
                   continue; // a cell is not the neighbor to itself
               }
               if (outOfBounds(lattice, x+xx, y+yy)) {
                   continue; // if it's out of bounds, don't check it
               }
              
               // only add it if there's no cell there
               if (lattice[x+xx][y+yy] == null) {
                   ArrayList<Integer> spot = new ArrayList<Integer>();
                   spot.add(x+xx);
                   spot.add(y+yy);
                  
                   freeSpots.add(spot);
               }
           }
       }
       return freeSpots; // return the array
   }
  
   /**
    * Return whether a position is out of bounds of a lattice.
    * @param lattice: a Cell[][] object that is the current state of the environment
    * @param x: the current x position.
    * @param y: the current y position.
    * @return boolean: true if it is out of bounds, false if it is not.
    */
   static boolean outOfBounds(Cell[][] lattice, int x, int y) {
       return (x < 0 || y < 0 || x >= lattice.length || y >= lattice[0].length);
   }
  
   /**
    * Prints out the current lattice in an easy to read structure.
    * @param lattice: a Cell[][] object that is the current state of the environment
    */
   static void print(Cell[][] lattice) {
       // for every single spot in the lattice
       for (int i = 0; i < lattice.length; i++) {
           for (int j = 0; j < lattice[0].length; j++) {
               if (lattice[i][j] == null) { // if it is null, print out O for empty spot
                   System.out.print("O ");
               } else if (!lattice[i][j].isStem) {
                   System.out.print("N, "); // print out N for NON-STEM cell
               } else {
                   System.out.print("S, "); // print out S for STEM cell
               }
           }
           System.out.println();
       }
   }
}


