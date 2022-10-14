package ch.irb.utils;

import java.io.File;
import java.text.Normalizer;
import java.util.*;
import java.util.Map.Entry;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
 * Copyright 2017 - Mathilde Foglierini Perez
 * This code is distributed open source under the terms of the GNU Free Documention License.
 * <p>
 * Utils class.
 */

public final class Consts {

    public static final String ls = System.getProperty("line.separator");
    public static final String fs = System.getProperty("file.separator");
    public static final String[] aminoAcidsList = {"A", "G", "K", "N", "S", "D", "V", "L", "T", "I", "E", "Q",
            "R", "H", "M", "P", "C", "F", "W", "Y", "X", "-"};
    @SuppressWarnings({"unchecked", "rawtypes"})
    public static final ArrayList<String> aminoAcids = new ArrayList(Arrays.asList(aminoAcidsList));

    public static final String[] immunoSeqVGeneList = {"TCRBV01-01", "TCRBV02-01", "TCRBV03", "TCRBV04", "TCRBV04-01", "TCRBV04-02"
            , "TCRBV04-03", "TCRBV04-04", "TCRBV05", "TCRBV05-01", "TCRBV05-02", "TCRBV05-03", "TCRBV05-04", "TCRBV05-05", "TCRBV05-06"
            , "TCRBV05-07", "TCRBV05-08", "TCRBV06", "TCRBV06-01", "TCRBV06-04", "TCRBV06-05", "TCRBV06-06", "TCRBV06-07", "TCRBV06-08"
            , "TCRBV06-09", "TCRBV07", "TCRBV07-01", "TCRBV07-02", "TCRBV07-03", "TCRBV07-04", "TCRBV07-05", "TCRBV07-06", "TCRBV07-07"
            , "TCRBV07-08", "TCRBV07-09", "TCRBV08-02", "TCRBV09-01", "TCRBV10-01", "TCRBV10-02", "TCRBV10-03", "TCRBV11", "TCRBV11-01"
            , "TCRBV11-02", "TCRBV11-03", "TCRBV12", "TCRBV12-01", "TCRBV12-02", "TCRBV12-05", "TCRBV13-01", "TCRBV14-01", "TCRBV15-01"
            , "TCRBV16-01", "TCRBV18-01", "TCRBV19-01", "TCRBV20", "TCRBV20-01", "TCRBV21-01", "TCRBV22-01", "TCRBV23-01", "TCRBV24"
            , "TCRBV25", "TCRBV25-01", "TCRBV27-01", "TCRBV28-01", "TCRBV29-01", "TCRBV30-01"};
    public static final ArrayList<String> immunoSeqVGenes = new ArrayList<String>(Arrays.asList(immunoSeqVGeneList));


    @SuppressWarnings("rawtypes")
    public static <K extends Comparable, V extends Comparable> Map<K, V> sortByValuesAsc(Map<K, V> map) {
        List<Map.Entry<K, V>> entries = new LinkedList<Map.Entry<K, V>>(map.entrySet());

        Collections.sort(entries, new Comparator<Map.Entry<K, V>>() {

            @SuppressWarnings("unchecked")
            public int compare(Entry<K, V> o1, Entry<K, V> o2) {
                return o1.getValue().compareTo(o2.getValue());
            }
        });

        // LinkedHashMap will keep the keys in the order they are inserted
        // which is currently sorted on natural ordering
        Map<K, V> sortedMap = new LinkedHashMap<K, V>();

        for (Map.Entry<K, V> entry : entries) {
            sortedMap.put(entry.getKey(), entry.getValue());
        }

        return sortedMap;
    }

    @SuppressWarnings("rawtypes")
    public static <K extends Comparable, V extends Comparable> Map<K, V> sortByValuesDesc(Map<K, V> map) {
        List<Map.Entry<K, V>> entries = new LinkedList<Map.Entry<K, V>>(map.entrySet());

        Collections.sort(entries, new Comparator<Map.Entry<K, V>>() {

            @SuppressWarnings("unchecked")
            public int compare(Entry<K, V> o1, Entry<K, V> o2) {
                return o2.getValue().compareTo(o1.getValue());
            }
        });

        // LinkedHashMap will keep the keys in the order they are inserted
        // which is currently sorted on natural ordering
        Map<K, V> sortedMap = new LinkedHashMap<K, V>();

        for (Map.Entry<K, V> entry : entries) {
            sortedMap.put(entry.getKey(), entry.getValue());
        }

        return sortedMap;
    }

    /*
     * Get the extension of a file.
     */
    public static String getExtension(File f) {
        String ext = null;
        String s = f.getName();
        int i = s.lastIndexOf('.');

        if (i > 0 && i < s.length() - 1) {
            ext = s.substring(i + 1).toLowerCase();
        }
        return ext;
    }

    /*
    To use this method the IDs have to be: genBank_ID|country|collection_date
     */
    public static HashMap<String, String> keepUniqueSequencesAndCreateId(HashMap<String, String> idToSeq_toCheck) {
        //Finally we keep only the unique sequence and change the ID accordingly
        HashMap<String, ArrayList<String>> seqToIds = new HashMap<String, ArrayList<String>>();
        for (String id : idToSeq_toCheck.keySet()) {
            String seq = idToSeq_toCheck.get(id);
            ArrayList<String> ids = new ArrayList<String>();
            if (seqToIds.containsKey(seq)) {
                ids = seqToIds.get(seq);
            }
            ids.add(id);
            seqToIds.put(seq, ids);
        }
        HashMap<String, String> idToSeq = new HashMap<String, String>();
        for (String seq : seqToIds.keySet()) {
            ArrayList<String> ids = seqToIds.get(seq);
            String newId = ids.get(0);
            if (ids.size() > 1) {
                TreeMap<Integer, String> yearToId = new TreeMap<Integer, String>();
                //System.out.println("SEQ dUPLI "+ids.size());
                for (String id : ids) {
                    String date = id.split("\\|")[2];
                    //System.out.println(date+" for id "+id);
                    int year = 3000; //for NA
                    String pattern = "(.*)(\\d{4})(.*)";
                    // Create a Pattern object
                    Pattern r = Pattern.compile(pattern);
                    // Now create matcher object.
                    Matcher m = r.matcher(date);
                    if (m.find()) { //we have the date
                        //System.out.println("Found value: " + m.group(2));
                        year = Integer.valueOf(m.group(2)).intValue();
                    }
                    yearToId.put(Integer.valueOf(year), id);
                }
                //then we keep the oldest year related id and we add the number of duplicates
                newId = yearToId.get(yearToId.firstKey()) + " (" + ids.size() + ")";
                System.out.println(newId);
            }
            idToSeq.put(newId, seq);
        }
        return idToSeq;
    }

    /*Calculate the geometric mean from a double array*/
    public static float geometricMean(double[] arr) {
        // declare sum variable and
        // initialize it to 1.
        float sum = 0;
        int n = arr.length;

        // Compute the sum of all the
        // elements in the array.
        for (int i = 0; i < n; i++)
            sum = sum + (float) Math.log(arr[i]);

        // compute geometric mean through formula
        // antilog(((log(1) + log(2) + . . . + log(n))/n)
        // and return the value to main function.
        sum = sum / n;

        return (float) Math.exp(sum);
    }


    public static double getStandardDeviation(double[] arr) {

        int n = arr.length;
        double standardDeviation = 0.0;
        double mean = 0.0;
        double res = 0.0;
        double sq = 0.0;
        double sum=0;

            for (int i = 0; i < n; i++) {
                sum = sum + arr[i];
            }
            mean = sum / (n);

        for (int i = 0; i < n; i++) {
           // standardDeviation = standardDeviation + Math.pow((arr[i] - mean), 2);
            standardDeviation = standardDeviation + (arr[i] - mean) * (arr[i] - mean);

        }
        sq = standardDeviation / n;
        res = Math.sqrt(sq);
        return res;
    }

    public static double getGeometricStandardDeviation(double[] arr, double mean) {

        int n = arr.length;
        double standardDeviation = 0.0;
        double res = 0.0;
        double res1 = 0.0;
        double sq = 0.0;

        for (int i = 0; i < n; i++) {
            double ratio = arr[i]/mean;
            standardDeviation = standardDeviation + Math.pow(Math.log(ratio),2);

        }

        sq = standardDeviation / n;
        res = Math.sqrt(sq);
        res1 = Math.exp(res);
        return res1;
    }

    public static String unaccent(String src) {
        return Normalizer
                .normalize(src, Normalizer.Form.NFD)
                .replaceAll("[^\\p{ASCII}]", "");
    }
}
