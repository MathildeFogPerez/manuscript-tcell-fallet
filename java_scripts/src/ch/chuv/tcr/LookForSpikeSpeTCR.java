package ch.chuv.tcr;

import ch.irb.utils.Consts;
import ch.irb.utils.IMGTnomenclatureConverter;

import java.io.*;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedHashMap;

/*
Copyright 2022 - Mathilde Foglierini Perez
This class will take in input the AIRR file produced by scirpy (cleaned TCR dataset after RNA and TCR processing)
It will then check in minigene, peptide class I and II files if there are some matches.
Those files come from MIRA release 2.1, all TCR sequences are Spike specific (total ~160'000 sequences).
It will also look for SARS-Cov-2 specific TCRs coming from vdjdb database.
There are 2 modes to search: CDR3 aa or CDR3 aa + V family. The latest was used for publication.
 */

public class LookForSpikeSpeTCR {

    private File inputAirrFile = new File("/MYPATH/manuscript-tcell-fallet/processed_files/10X_airrfile_data_cleaned.tsv");
    //here set the path to the two use databases: MIRA and vdjdb
    private File vdjdbFile = new File("/MYPATH/vdjdb/vdjdb_sarscov_27_10_21.tsv");
    private File miradbFiles = new File("/MYPATH/ImmuneCODE-MIRA-Release002.1/");
    private String mode ="cdr3_Vfamily"; //can be cdr3 or cdr3_Vfamily if we include the v FAMILY
    private LinkedHashMap<String,TCR> sequenceIdToTcr = new LinkedHashMap<>();
    private HashMap<String,Integer> sequenceIdToHits = new HashMap<>();
    private HashMap<String,HashSet<String>> sequenceIdToExps = new HashMap<>();
    private HashMap<String,HashSet<String>> sequenceIdToORF = new HashMap<>();
    private String subjectHeader = null;
    private HashMap<String,String> expToLine = new HashMap<>();
    private int vdjdbHits=0;

    public static void main(String[] args) throws IOException {
        new LookForSpikeSpeTCR();
    }

    public LookForSpikeSpeTCR() throws IOException {
        parseData();
        parseSubjectMetadata();
        for (File file: miradbFiles.listFiles()){
            if (file.getName().contains("detail")){
                findMiraMatches(file);
            }
        }
        findVdjdbMatches();
        writeFinalOut();
    }

    private void parseData() throws IOException {
        HashMap<String,Integer> headerToIndex = new HashMap<>();
        BufferedReader fileReader = new BufferedReader(new FileReader(inputAirrFile));
        String line="";
        int index=0;
        while ((line= fileReader.readLine())!=null){
            String[] cells =line.split("\t");
            if (index==0){
                int i=0;
                for (String cell: cells){
                    headerToIndex.put(cell,i);
                    i++;
                }
            }
            else{
                String cellbarcode = cells[headerToIndex.get("cell_id")];
                String sequence_id= cells[headerToIndex.get("sequence_id")];
                TCR tcr;
                if (sequenceIdToTcr.containsKey(sequence_id)){
                    tcr= sequenceIdToTcr.get(sequence_id);
                }
                else {
                    tcr = new TCR(cellbarcode,sequence_id);
                }
                tcr.setDonor(cells[headerToIndex.get("sample")]);
                tcr.setLeiden(cells[headerToIndex.get("leiden")]);
                tcr.setClone_id(cells[headerToIndex.get("clone_id")]);
                tcr.setClone_id_size(cells[headerToIndex.get("clone_id_size")]);
                if (cells[headerToIndex.get("locus")].equals("TRA")){
                    tcr.setCdr3a(cells[headerToIndex.get("junction_aa")]);
                    tcr.setVa(cells[headerToIndex.get("v_call")]);
                    tcr.setJa(cells[headerToIndex.get("j_call")]);
                }else{
                    tcr.setCdr3b(cells[headerToIndex.get("junction_aa")]);
                    tcr.setVb(cells[headerToIndex.get("v_call")]);
                    tcr.setJb(cells[headerToIndex.get("j_call")]);
                }
                sequenceIdToTcr.put(sequence_id,tcr);
            }
            index++;
        }
        fileReader.close();
        System.out.println(" Number of UNPAIRED TCRs processed from the AIRR file: "+sequenceIdToTcr.size());
    }

    private void parseSubjectMetadata() throws IOException  {
        File file = new File(miradbFiles.getPath()+ Consts.fs+"subject-metadata.csv");
        BufferedReader fileReader = new BufferedReader(new FileReader(file));
        String line="";
        int index=0;
        while ((line= fileReader.readLine())!=null){
            String lineToProcess= line.replaceAll(",","\t");
            if (index==0){
                subjectHeader =lineToProcess;
            }
            else{
                String exp = lineToProcess.split("\t")[0].replaceAll("\"","");
                expToLine.put(exp,lineToProcess.replaceAll("\"",""));
            }
            index++;
        }
        fileReader.close();
        System.out.println("Number of subject in MIRA database: "+expToLine.size());
    }

    private void findMiraMatches(File miraFile) throws IOException  {
        System.out.println("Proces Mira file "+miraFile.getName());
        HashMap<String,Integer> headerToIndex = new HashMap<>();
        File outfile = new File(inputAirrFile.getParent()+Consts.fs
                + inputAirrFile.getName().replace(".tsv","_")
                +miraFile.getName().replace(".csv","")
                +"_"+mode+".tsv");
        BufferedWriter out = new BufferedWriter(new FileWriter(outfile));
        BufferedReader fileReader = new BufferedReader(new FileReader(miraFile));
        String line="";
        String header=null;
        int index=0;
        while ((line= fileReader.readLine())!=null){
            String lineToProcess= line.replaceAll(",","\t").replaceAll("\"","");
            String[] cells =lineToProcess.split("\t");
            if (index==0){
                header=lineToProcess;
                out.write("sequence_id\tdonor\tleiden\tbarcode\tvb\tjb\tcdr3b\t"+header+"\t"+subjectHeader+Consts.ls);
                int i=0;
                for (String cell: cells){
                    String cellToStore =cell.replaceAll("\"","");
                    headerToIndex.put(cellToStore,i);
                    //System.out.println(cellToStore);
                    i++;
                }
            }
            else{
                //in the case of CDR3+Vgene, DO NOT FORGET to convert to Immunoseq nomenclature
                String bioIdent= cells[headerToIndex.get("TCR BioIdentity")];
                String cdr3=bioIdent.split("\\+")[0];
                String v=bioIdent.split("\\+")[1].split("-")[0]; //we take the V family, EASIER to check
                String toCheck1= cdr3;
                if (!mode.equals("cdr3")){
                    toCheck1+=v;
                }
                //go through each clonotypes
                for (String seqId: sequenceIdToTcr.keySet()) {
                    //Do only TRB!!!!
                    if (seqId.contains("_TRB_")){
                        TCR tcr = sequenceIdToTcr.get(seqId);
                        String toCheck2 = tcr.getCdr3b();
                        if (!mode.equals("cdr3")) {
                            String v2 = new IMGTnomenclatureConverter(tcr.getVb(), true).getImmunoSeqFamily();
                            toCheck2 += v2;
                        }
                        //System.out.println("CHECK MIRA "+toCheck1+ " with "+toCheck2);
                        if (toCheck1.equals(toCheck2)) {
                            System.out.println("--> HIT Mira " + toCheck1 + " with " + toCheck2);
                            int count =1;
                            if (sequenceIdToHits.containsKey(seqId)){
                                count += sequenceIdToHits.get(seqId);
                            }
                            sequenceIdToHits.put(seqId,count);
                            String exp = cells[headerToIndex.get("Experiment")];
                            HashSet<String> exps = new HashSet<>();
                            if (sequenceIdToExps.containsKey(seqId)){
                                exps= sequenceIdToExps.get(seqId);
                            }
                            exps.add(exp);
                            sequenceIdToExps.put(seqId,exps);
                            String orf = cells[headerToIndex.get("ORF Coverage")];
                            HashSet<String> orfs = new HashSet<>();
                            if (sequenceIdToORF.containsKey(seqId)){
                                orfs= sequenceIdToORF.get(seqId);
                            }
                            orfs.add(orf);
                            sequenceIdToORF.put(seqId, orfs);
                            out.write(tcr.getSequence_id() + "\t" + tcr.getDonor() + "\t" + tcr.getLeiden() + "\t" + tcr.getCellbarcode() + "\t" + tcr.getVb()
                                + "\t" + tcr.getJb() + "\t" + tcr.getCdr3b() + "\t" + lineToProcess + "\t" + expToLine.get(exp) + Consts.ls);
                        }
                    }
                }
            }
            index++;
        }
        fileReader.close();
        out.close();
        System.out.println();
    }

    private void findVdjdbMatches() throws IOException {
        System.out.println("Proces VDJDB file "+vdjdbFile.getName());
        HashMap<String,Integer> headerToIndex = new HashMap<>();
        File outfile = new File(inputAirrFile.getParent()+Consts.fs
                + inputAirrFile.getName().replace(".tsv","_")
                +vdjdbFile.getName().replace(".tsv","")
                +"_"+mode+".tsv");
        BufferedWriter out = new BufferedWriter(new FileWriter(outfile));
        BufferedReader fileReader = new BufferedReader(new FileReader(vdjdbFile));
        String line="";
        int index=0;
        String header=null;
        while ((line= fileReader.readLine())!=null){
            String[] cells =line.replaceAll("\"","").split("\t");
            if (index==0){
                header=line;
                out.write("sequence_id\tdonor\tleiden\tbarcode\tva/b\tja/b\tcdr3a/b\t"+header+"\t"+subjectHeader+Consts.ls);
                int i=0;
                for (String cell: cells){
                    String cellToStore =cell.replaceAll("\"","");
                    headerToIndex.put(cellToStore,i);
                    //System.out.println(cellToStore);
                    i++;
                }
            }
            else{
                String cdr3=cells[headerToIndex.get("CDR3")];
                String v=getVFamily(cells[headerToIndex.get("V")]); //we take the V family, EASIER to check
                String tr=cells[headerToIndex.get("Gene")];
                String toCheck1= cdr3;
                if (!mode.equals("cdr3")){
                    toCheck1+=v;
                }
                for (String seqId: sequenceIdToTcr.keySet()) {
                    //Do only TRA or TRB
                    if (seqId.contains("_"+tr+"_")){
                        TCR tcr = sequenceIdToTcr.get(seqId);
                        String toCheck2;
                        if (tr.equals("TRA")){
                            toCheck2 = tcr.getCdr3a();
                        }
                        else{
                            toCheck2 = tcr.getCdr3b();
                        }
                        if (!mode.equals("cdr3")) {
                            String v2;
                            if (tr.equals("TRA")){
                                v2 = getVFamily(tcr.getVa());//get family only
                            }else{
                                v2 = getVFamily(tcr.getVb());//get family only
                            }
                            toCheck2 += v2;
                        }
                       // System.out.println("CHECK VDJDB "+toCheck1+ " with "+toCheck2);
                        if (toCheck1.equals(toCheck2)) {
                            System.out.println("--> HIT VDJDB " + toCheck1 + " with " + toCheck2);
                            int count =1;
                            if (sequenceIdToHits.containsKey(seqId)){
                                count += sequenceIdToHits.get(seqId);
                            }
                            sequenceIdToHits.put(seqId,count);
                            vdjdbHits++;
                            HashSet<String> exps = new HashSet<>();
                            if (sequenceIdToExps.containsKey(seqId)){
                                exps = sequenceIdToExps.get(seqId);
                            }
                            //get donor if exists
                            String meta= cells[headerToIndex.get("Meta")];
                            String[] metas = meta.replaceAll("\"","").split(",");
                            String donor="NA";
                            for (String met: metas){
                                if (met.contains("subject.id")){
                                    donor=met.replace("subject.id:","");
                                    if (donor.length()==0){
                                        donor="NA";
                                    }
                                }
                            }
                            exps.add(donor);
                            sequenceIdToExps.put(seqId,exps);
                            out.write(tcr.getSequence_id() + "\t" + tcr.getDonor() + "\t" + tcr.getLeiden() + "\t" + tcr.getCellbarcode() + "\t" );
                            if (tr.equals("TRB")) {
                                out.write(tcr.getVb() + "\t" + tcr.getJb() + "\t" + tcr.getCdr3b());
                            }
                            else{
                                out.write(tcr.getVa() + "\t" + tcr.getJa() + "\t" + tcr.getCdr3a());
                            }
                            out.write("\t" + line + Consts.ls);
                        }
                    }
                }
            }
            index++;
        }
        fileReader.close();
        out.close();
    }

    private String getVFamily (String vGene){
        String vFamily = vGene;
        if (vGene.contains("-")){
            vFamily= vGene.split("-")[0];
        }
        else if (vGene.contains("*")){
            vFamily = vGene.split("\\*")[0];
        }
        return vFamily;
    }

    private void writeFinalOut() throws IOException {
        BufferedWriter out2 = new BufferedWriter(new FileWriter(new File(inputAirrFile.getParent()+Consts.fs+"cell_barcode_to_sarcov2_spe.tsv")));
        BufferedWriter out = new BufferedWriter(new FileWriter(new File(inputAirrFile.getParent()+Consts.fs
                +inputAirrFile.getName().replace(".tsv","_"+mode+"_withHITs_v3.tsv"))));
        out.write("sequence_id\tdonor\tleiden\tclone_id\tclone_id_size\tbarcode\tva/b\tja/b\tcdr3a/b\tisSarsCov2Spe" +
                "\thits_number\texps_number\texps\tSARSCOV2.ORF"+Consts.ls);
        for (String seqId: sequenceIdToTcr.keySet()){
            TCR tcr = sequenceIdToTcr.get(seqId);
            if (seqId.contains("_TRB_")) {
                out.write(tcr.getSequence_id() + "\t" + tcr.getDonor() + "\t" + tcr.getLeiden() + "\t" +tcr.getClone_id()+ "\t" +
                        tcr.getClone_id_size()+"\t" + tcr.getCellbarcode() + "\t" + tcr.getVb()
                        + "\t" + tcr.getJb() + "\t" + tcr.getCdr3b() + "\t" );
                if (sequenceIdToHits.containsKey(seqId)){
                    out.write("sarscov2\t"+sequenceIdToHits.get(seqId)+"\t"+sequenceIdToExps.get(seqId).size()+"\t"+sequenceIdToExps.get(seqId)
                            + "\t" +sequenceIdToORF.get(seqId));
                    out2.write(tcr.getCellbarcode() + "\t" +sequenceIdToORF.get(seqId)+Consts.ls);
                }else{
                    out.write("unknown\t0\t0\t0");
                }
                out.write(Consts.ls);
            }
        }
        out.close();
        out2.close();
        System.out.println("\nIn mode "+mode+",  Total number of TCR with Spike Spe: "+sequenceIdToHits.size());
        System.out.println("Number of TCR VDJDB hits: "+vdjdbHits);
        //int count = sequenceIdToHits.size()-vdjdbHits;
        //System.out.println("Number of TCR Spike specific to consider further: "+count);
    }



    public class TCR{
        private String sequence_id;
        private String cellbarcode;
        private String clone_id;
        private String clone_id_size;
        private String donor;
        private String va;
        private String vb;
        private String ja;
        private String jb;
        private String cdr3a;
        private String cdr3b;
        private String leiden;
        private String specificity="NA"; //can be spike
        private int tra_hits=0;
        private int trb_hits=0;

        public TCR(String cellbarcode,String sequence_id){
            this.cellbarcode=cellbarcode;
            this.sequence_id = sequence_id;
        }

        public String getCellbarcode() {
            return cellbarcode;
        }

        public void setCellbarcode(String cellbarcode) {
            this.cellbarcode = cellbarcode;
        }

        public String getClone_id() {
            return clone_id;
        }

        public void setClone_id(String clone_id) {
            this.clone_id = clone_id;
        }

        public String getClone_id_size() {
            return clone_id_size;
        }

        public void setClone_id_size(String clone_id_size) {
            this.clone_id_size = clone_id_size;
        }

        public String getSequence_id() {
            return sequence_id;
        }

        public String getDonor() {
            return donor;
        }

        public void setDonor(String donor) {
            this.donor = donor;
        }

        public String getVa() {
            return va;
        }

        public void setVa(String va) {
            this.va = va;
        }

        public String getVb() {
            return vb;
        }

        public void setVb(String vb) {
            this.vb = vb;
        }

        public String getJa() {
            return ja;
        }

        public void setJa(String ja) {
            this.ja = ja;
        }

        public String getJb() {
            return jb;
        }

        public void setJb(String jb) {
            this.jb = jb;
        }

        public String getCdr3a() {
            return cdr3a;
        }

        public void setCdr3a(String cdr3a) {
            this.cdr3a = cdr3a;
        }

        public String getCdr3b() {
            return cdr3b;
        }

        public void setCdr3b(String cdr3b) {
            this.cdr3b = cdr3b;
        }

        public String getLeiden() {
            return leiden;
        }

        public void setLeiden(String leiden) {
            this.leiden = leiden;
        }

        public String getSpecificity() {
            return specificity;
        }

        public void setSpecificity(String specificity) {
            this.specificity = specificity;
        }

        public int getTra_hits() {
            return tra_hits;
        }

        public void setTra_hits(int tra_hits) {
            this.tra_hits = tra_hits;
        }

        public int getTrb_hits() {
            return trb_hits;
        }

        public void setTrb_hits(int trb_hits) {
            this.trb_hits = trb_hits;
        }


    }

}
