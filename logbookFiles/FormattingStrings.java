package logbookFiles;

import java.util.ArrayList;
import java.util.LinkedHashSet;
import java.util.Scanner;
import java.util.Set;

public class FormattingStrings {
   public static void main(String[] args) {
       Scanner in = new Scanner(System.in);
       ArrayList<String> lines = new ArrayList<>();
      
       while (in.hasNext()) {
           lines.add(in.nextLine());
           if (in.next().equals("qqq")) {
               break;
           }
       }
       ArrayList<String> genes = new ArrayList<>();
      
       for (int i = 0; i < lines.size(); i++) {
           String entry = lines.get(i);
           int start = -1;
           int end = -1;
          
           for (int j = 0; j < entry.length(); j++) {
               if (entry.charAt(j) == '"' && start == -1) {
                   start = j;
               } else if (entry.charAt(j) == '"') {
                   end = j;
                   genes.add(entry.substring(start+1, end));
                   start = -1;
               }
           }
       }  
       ArrayList<String> unique = removeDuplicates(genes);
       for (int i = 0; i < unique.size(); i++) {
           System.out.print(unique.get(i) + "\n");
       }
       in.close();
   }
  
   // Function to remove duplicates from an ArrayList
   public static <T> ArrayList<T> removeDuplicates(ArrayList<T> list)
   {
        // Create a new LinkedHashSet
       Set<T> set = new LinkedHashSet<>();
        // Add the elements to set
       set.addAll(list);
        // Clear the list
       list.clear();
        // add the elements of set
       // with no duplicates to the list
       list.addAll(set);
        // return the list
       return list;
   }
}
