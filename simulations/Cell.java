package simulations;

import java.util.Random;
 
public class Cell {
   int x, y; // position of the cell in the grid
   int timeToNextDivision; // amount of time until next division
   boolean isStem;
   double symmetricDivision; // probability that a cell divides symmetrically
   int telomereLength; // length of a telomere that indicates when a cell dies
   double spontaneousDeath; // probability that a cell dies at random
  
   // initializing all values
   Cell(int x, int y, boolean stem, double tLength, double division, double death) {
       this.x = x;
       this.y = y;
       Random rand = new Random();
       timeToNextDivision = (int) rand.nextDouble() * 60;
       isStem = stem;
       symmetricDivision = division;
       telomereLength = 10;
       spontaneousDeath = death;
   }
  
   // print out the Cell with its type and position
   public String toString() {
       if (isStem) return "S " + x + " " + y;
       return "N " + x + " " + y;
   }
}


