/** **************************************
 * Name:		Tyler Haddox
 * Username:	NA
 * Problem Set: PS3
 * Due Date: March 3, 2020
 ************************************** */
package aips4;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.BitSet;
import java.util.Comparator;
import java.util.Random;
import java.util.TreeSet;

public class AI_ps4_Tools {

    private TreeSet<Trait> popSel = new TreeSet<>(new PopulationComparator());
    private ArrayList<Gene> genes = new ArrayList<>();
    private Trait savedMaxTrait = new Trait();
    private Random r = new Random();
    private int sampleSize = 0;
    private int popCnt = 11;
    private int maxItrCnt = 10;
    private float ksCapacity = 100;
    private float mutPer = 0.01f;
    private float minFitScore = 300;
    private int localPeakCnt = 0;
  
        
    public ArrayList<String> improvedGenetic() {
        ArrayList<String> rstr = new ArrayList<>();
        savedMaxTrait = new Trait();
        int diff = getMaxItrCnt() + 1;
        while (getMaxItrCnt() > 0 && getMinFitScore() > popSel.last().v) {            
            Trait p1 = getImprovedRandomTrait(2);
            Trait p2 = getRandomTrait();
            breed(p1, p2, 3);
            removeLow(3);
            rstr.add(improvedMeanVariance(diff));             
            setMaxItrCnt(getMaxItrCnt() - 1);
        }
        if(popSel.last().v > savedMaxTrait.v){
            savedMaxTrait = popSel.last();
        }
        rstr.add(printMaxFit(savedMaxTrait));        
//        try {
//            writeToFile(savedMaxTrait.v);
//        } catch (IOException ex) {
//            Logger.getLogger(AI_ps4_Tools.class.getName()).log(Level.SEVERE, null, ex);
//        }              
        return rstr;
    }
    
        public ArrayList<String> genetic() {
        ArrayList<String> rstr = new ArrayList<>();
        int diff = getMaxItrCnt() + 1;
        while (getMaxItrCnt() > 0 && getMinFitScore() > popSel.last().v) {
            Trait p1 = getRandomTrait();
            Trait p2 = getRandomTrait();
            breed(p1, p2, 3);
            removeLow(3);
            rstr.add(meanVariance(diff));  
            setMaxItrCnt(getMaxItrCnt() - 1);
        }
        rstr.add(printMaxFit(popSel.last()));        
//        try {
//            writeToFile(popSel.last().v);
//        } catch (IOException ex) {
//            Logger.getLogger(AI_ps4_Tools.class.getName()).log(Level.SEVERE, null, ex);
//        }        
        return rstr;
    }
        
    public String improvedMeanVariance(int diff){
        float mean = calcMean();
        float var = calcVar(mean);
        int itr = diff - getMaxItrCnt();
        String mv = String.format("%d:\tmean=%.0f    var=%.0f\n", itr, mean, var);
        System.out.printf(mv);  
        if(var < 1){
            localPeakCnt++;
        }else{
            localPeakCnt = 0;
        }
        if(localPeakCnt == 3){
            if(savedMaxTrait.v < popSel.last().v){
                savedMaxTrait = popSel.last();
            }
            popSel.clear();
            initialSelection();
        }        
        return mv;    
    }
    
    public String meanVariance(int diff){
        float mean = calcMean();
        float var = calcVar(mean);
        int itr = diff - getMaxItrCnt();
        String mv = String.format("%d:\tmean=%.0f    var=%.0f\n", itr, mean, var);
        System.out.printf(mv);        
        return mv;    
    }
    
    public float calcVar(float mean){
        float var = 0;
        for (Trait trait : popSel) {
            var += Math.pow(trait.v - mean, 2);
        }
        var = var / (getPopCnt() - 1);        
        return var;
    }
    
    public float calcMean(){
        float mean = 0;
        for (Trait trait : popSel) {
            mean += trait.v;
        }
        mean = mean / getPopCnt();
        return mean;
    }
        
    public Trait getImprovedRandomTrait(int k) {        
        Object[] obj = popSel.toArray();
        int high = r.nextInt(obj.length);
        if(high < obj.length/k){
            high += obj.length/k;
        }        
        return (Trait) obj[high];        
    }
        
    public String printMaxFit(Trait maxFitTrait) {
        String traitBits = "";
        for (int i = getSampleSize() - 1; i > -1; i--) {
            if (maxFitTrait.b.get(i)) {
                traitBits = traitBits.concat("1");
            } else {
                traitBits = traitBits.concat("0");
            }
        }        
        System.out.printf("\nTotal Weight:\t%.0f\nTotal Value:\t$ %.2f\n"
                + "Layout:\t\t%s\n", maxFitTrait.w, maxFitTrait.v, traitBits);
        
        return String.format("\nTotal Weight:\t%.0f\nTotal Value:\t$ %.2f\nLayout:\t%s", maxFitTrait.w, maxFitTrait.v, traitBits);
    }    
        
    public Trait getRandomTrait() {
        Object[] obj = popSel.toArray();
        return (Trait) obj[r.nextInt(obj.length)];
    }
       
    public void removeLow(int i) {
        while (i > 0) {
            popSel.pollFirst();
            i--;
        }
    }

    public void breed(Trait p1, Trait p2, int k) {
        BitSet b = new BitSet(getSampleSize());
        while (k > 0) {
            Trait trait = new Trait();
            BitSet bChild = new BitSet(getSampleSize());
            do {
                b.clear();
                bChild.clear();
                b = setRandomBits(b, 1);
                for (int i = 0; i < getSampleSize(); i++) {
                    if (b.get(i)) {
                        if (p1.b.get(i)) {
                            bChild.set(i);
                        }
                    } else {
                        if (p2.b.get(i)) {
                            bChild.set(i);
                        }
                    }
                }
                bChild = mutate(bChild);
                trait = new Trait(bChild);
            } while (trait.w > getKsCapacity());
            popSel.add(trait);
            k--;
        }
    }

    public void improvedInitialSelection() {
        int i = 0;
        while (i < getPopCnt()) {
            Trait trait = new Trait();
            BitSet b = new BitSet(getSampleSize());
            do {
                b.clear();
                b = improvedSetRandomBits(b, 2);
                trait = new Trait(b);
            } while (trait.w > getKsCapacity());
            popSel.add(trait);
            i++;
        }
    }    

    public void initialSelection() {
        int i = 0;
        while (i < getPopCnt()) {
            Trait trait = new Trait();
            BitSet b = new BitSet(getSampleSize());
            do {
                b.clear();
                b = setRandomBits(b, 2);
                trait = new Trait(b);
            } while (trait.w > getKsCapacity());
            popSel.add(trait);
            i++;
        }
    }
    
    public BitSet mutate(BitSet bChild) {
        int m = (int) (getSampleSize() * getMutPer());
        if(m < 1){
            m = 1;
        }
        while (m > 0) {            
            bChild.flip(r.nextInt(getSampleSize()));
            m--;
        }
        return bChild;
    }

    public BitSet improvedSetRandomBits(BitSet b, int size) {
        int bits = getSampleSize() / size;
        for (int i = 0; i < bits; i++) {
            int ran = r.nextInt(getSampleSize());
            if(ran < bits){
                ran = ran + bits;
            }
            b.set(ran);
        }
        return b;
    }
    
    public BitSet setRandomBits(BitSet b, int size) {
        for (int i = 0; i < getSampleSize() / size; i++) {
            b.set(r.nextInt(getSampleSize()));
        }
        return b;
    }

    public void readFile(String fn) {
        File f = new File(fn);
        String strArray[];
        String read;
        try {
            BufferedReader br = new BufferedReader(new FileReader(f));
            while ((read = br.readLine()) != null) {
                strArray = read.split(",");
                Gene item = new Gene(Float.parseFloat(strArray[0]), Float.parseFloat(strArray[1]));
                genes.add(item);
            }
            br.close();
        } catch (FileNotFoundException e) {
            System.out.println("file not found");
        } catch (NumberFormatException | IOException e) {
        }
        genes.trimToSize();
        genes.sort(new MyComparator());
        setSampleSize(genes.size());
    }
    
    public void writeToFile(float v) throws IOException{
        BufferedWriter out = new BufferedWriter(
                new FileWriter("testResults.txt", true));
        String str = String.format("%.2f\n", v); 
        out.write(str);
        out.close();
    }    
    
    public int getSampleSize() {
        return sampleSize;
    }

    public void setSampleSize(int sampleSize) {
        this.sampleSize = sampleSize;
    }

    public int getPopCnt() {
        return popCnt;
    }

    public void setPopCnt(int popCnt) {
        this.popCnt = popCnt;
    }

    public int getMaxItrCnt() {
        return maxItrCnt;
    }

    public void setMaxItrCnt(int maxItrCnt) {
        this.maxItrCnt = maxItrCnt;
    }

    public float getKsCapacity() {
        return ksCapacity;
    }

    public void setKsCapacity(float ksCapacity) {
        this.ksCapacity = ksCapacity;
    }

    public float getMutPer() {
        return mutPer;
    }

    public void setMutPer(float mutPer) {
        this.mutPer = mutPer/100;
    }

    public float getMinFitScore() {
        return minFitScore;
    }

    public void setMinFitScore(float minFitScore) {
        this.minFitScore = minFitScore;
    }

    public void resetAll(){
        genes.clear();
        popSel.clear();        
    }
        
    class MyComparator implements Comparator<Gene> {
        public int compare(Gene a, Gene b) {
            if (a.localFitVal > b.localFitVal) {
                return 1;
            } else if (a.localFitVal < b.localFitVal) {
                return -1;
            } else if (a.w > b.w) {
                return 1;
            } else if (a.w < b.w) {
                return -1;
            }
            return -1;
        }
    }

    class PopulationComparator implements Comparator<Trait> {
        public int compare(Trait a, Trait b) {
            if (a.v > b.v) {
                return 1;
            } else if (a.v < b.v) {
                return -1;
            } else if (a.w > b.w) {
                return 1;
            } else if (a.w < b.w) {
                return -1;
            }
            return -1;
        }
    }

    public class Trait extends Gene {
        BitSet b = new BitSet();
        
        public Trait() {}
        public Trait(float w, float v, BitSet b) {
            this.w = w;
            this.v = v;
            this.b = b;
        }
        public Trait(BitSet b) {
            this.b = b;
            for (int i = b.nextSetBit(0); i >= 0; i = b.nextSetBit(i + 1)) {
                if (i == Integer.MAX_VALUE) {
                    break;
                }
                Gene gene = genes.get(i);
                this.w += gene.w;
                this.v += gene.v;
            }
        }
    }

    public class Gene {
        float w = 0;
        float v = 0;
        float localFitVal;

        public Gene() {}
        public Gene(float w, float v) {
            this.w = w;
            this.v = v;
            this.localFitVal = v/w;
        }
    }

}
