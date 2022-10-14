package ch.chuv.tcr;

import ch.irb.utils.Consts;

import java.io.*;
import java.util.*;

/*
This class will take in input the AIRR file produced by scirpy (cleaned TCR dataset after RNA and TCR processing)
It will check in BD dataset if we find shared clonotypes
There are 2 modes to search: CDR3 aa or CDR3 aa + V family. The latest was used for publication.
 */

public class LookForTCRintoBDdataset {

    private File inputAirrFile = new File("/MYPATH/manuscript-tcell-fallet/analysis_out_10X/files/airrfile_data_cleaned.tsv");

    private File bdFile = new File("/MYPATH/manuscript-tcell-fallet/analysis_out_BD/files/airrfile_data_cleaned.tsv");
    private File sampleTagFile = new File("/MYPATH/manuscript-tcell-fallet/BD_WTA_pipeline_out/" +
            "BD_Sample_Tag_Calls.csv");

    private File cellTypeFile = new File("/MYPATH/manuscript-tcell-fallet/BD_WTA_pipeline_out/"
            + "BD_cell_type_experimental.csv");

    private File sarcov2File = new File("/MYPATH/manuscript-tcell-fallet/processed_files/"
            + "10X_cell_barcode_to_sarcov2_spe.tsv");
    private String mode = "cdr3_Vfamily"; //can be cdr3 or cdr3_Vfamily if we include the v FAMILY
    private LinkedHashMap<String, TCR> barcodeIdToTcr10X ;
    private LinkedHashMap<String, TCR> barcodeIdToTcrBD;
    private HashMap<String, Integer> sequenceIdToHits = new HashMap<>();
    private HashMap<String, Integer> sequenceIdToHitsWithTRA = new HashMap<>();
    private HashMap<String, String> cellIndexToSample = new HashMap<>();
    private HashMap<String, String> cellIndexToCellType = new HashMap<>();
    private HashMap<String, String> cellIndexToSarscov2 = new HashMap<>();


    public static void main(String[] args) throws IOException {
        new LookForTCRintoBDdataset();
    }

    public LookForTCRintoBDdataset() throws IOException {
        //parse
        barcodeIdToTcr10X = parseData(inputAirrFile);
        barcodeIdToTcrBD = parseData(bdFile);
        parseSampleTags();
        parseCellType();
        parseSarscov2matches();
        findBDmatches();
    }

    private LinkedHashMap<String, TCR>  parseData(File inputAirrFile ) throws IOException {
        LinkedHashMap<String, TCR> sequenceIdToTcr = new LinkedHashMap<>();
        HashMap<String, Integer> headerToIndex = new HashMap<>();
        BufferedReader fileReader = new BufferedReader(new FileReader(inputAirrFile));
        String line = "";
        int index = 0;
        while ((line = fileReader.readLine()) != null) {
            String[] cells = line.split("\t");
            if (index == 0) {
                int i = 0;
                for (String cell : cells) {
                    headerToIndex.put(cell, i);
                    i++;
                }
            } else {
                String cellbarcode = cells[headerToIndex.get("cell_id")];
                String sequence_id = cells[headerToIndex.get("sequence_id")];
                TCR tcr;
                if (sequenceIdToTcr.containsKey(cellbarcode)) {
                    tcr = sequenceIdToTcr.get(cellbarcode);
                } else {
                    tcr = new TCR(cellbarcode);
                }
                tcr.setDonor(cells[headerToIndex.get("sample")]);
                //tcr.setLeiden(cells[headerToIndex.get("leiden")]);
                tcr.setCell_type(cells[headerToIndex.get("cell_type")]);
                tcr.setClone_id(cells[headerToIndex.get("clone_id")]);
                tcr.setClone_id_size(cells[headerToIndex.get("clone_id_size")]);
                if (cells[headerToIndex.get("locus")].equals("TRA")) {
                    tcr.setSequence_id_A(sequence_id);
                    tcr.setCdr3a(cells[headerToIndex.get("junction_aa")]);
                    tcr.setVa(cells[headerToIndex.get("v_call")]);
                    tcr.setJa(cells[headerToIndex.get("j_call")]);
                    tcr.setVdj_a_aa_seq(cells[headerToIndex.get("sequence_aa")]);
                } else {
                    tcr.setSequence_id_B(sequence_id);
                    tcr.setCdr3b(cells[headerToIndex.get("junction_aa")]);
                    tcr.setVb(cells[headerToIndex.get("v_call")]);
                    tcr.setJb(cells[headerToIndex.get("j_call")]);
                    tcr.setVdj_b_aa_seq(cells[headerToIndex.get("sequence_aa")]);
                }
                sequenceIdToTcr.put(cellbarcode, tcr);
            }
            index++;
        }
        fileReader.close();
        System.out.println(" Number of UNPAIRED TCRs processed from the AIRR file: " + sequenceIdToTcr.size());
        return (sequenceIdToTcr);
    }

    private void parseSampleTags() throws IOException {
        BufferedReader fileReader = new BufferedReader(new FileReader(sampleTagFile));
        String line = "";
        while ((line = fileReader.readLine()) != null) {
            String[] cells = line.split(",");
            if (!line.matches("#.*") && !line.matches("Cell_Index.*")) {
                //System.out.println(line);
                cellIndexToSample.put(cells[0], cells[2]);
            }
        }
        fileReader.close();
    }

    private void parseCellType() throws IOException {
        BufferedReader fileReader = new BufferedReader(new FileReader(cellTypeFile));
        String line = "";
        while ((line = fileReader.readLine()) != null) {
            String[] cells = line.split(",");
            if (!line.matches("#.*") && !line.matches("Cell_Index.*")) {
                //System.out.println(line);
                cellIndexToCellType.put(cells[0], cells[1]);
            }
        }
        fileReader.close();
    }

    private void parseSarscov2matches() throws IOException {
        BufferedReader fileReader = new BufferedReader(new FileReader(sarcov2File));
        String line = "";
        while ((line = fileReader.readLine()) != null) {
            String[] cells = line.split("\t");
            if (!line.matches("#.*") && !line.matches("Cell_Index.*")) {
                //System.out.println(line);
                cellIndexToSarscov2.put(cells[0], cells[1]);
            }
        }
        fileReader.close();
    }

    private void findBDmatches() throws IOException {
        HashSet<String> bdCells = new HashSet<>();
        HashSet<String> tenxCells = new HashSet<>();
        HashSet<String> lineForSankeyGraphBP01 = new HashSet<>();
        HashSet<String> lineForSankeyGraphBP02 = new HashSet<>();
        System.out.println("Process DB file " + bdFile.getName());
        HashMap<String, Integer> headerToIndex = new HashMap<>();
        File outfile = new File(inputAirrFile.getParent() + Consts.fs
                + inputAirrFile.getName().replace(".tsv", "_BDmatches_FINAL.tsv"));
        BufferedWriter out = new BufferedWriter(new FileWriter(outfile));
        BufferedReader fileReader = new BufferedReader(new FileReader(bdFile));
        int tcrAhits=0;
        int missingTcrA=0;
        int noTcrAhits=0;
        String line = "";
        int index = 0;
        while ((line = fileReader.readLine()) != null) {
            String[] cells = line.split("\t");
            if (index == 0) {
                out.write("sequence_id_10X\tclone_id_10X\tclone_id_size_10X\tdonor_10X\tcell_type_10X\tbarcode_10X\tvb_10X\tjb_10X\tcdr3b_10X" +
                        "\tva_10X\tja_10X\tcdr3a_10X\tsarscov2.spe\tsample_tag_BD\tcell_type_BD\tva_BD\tja_BD\tcdr3a_BD");
                //add _BD to headers
                for (String cell:cells){
                    out.write("\t"+cell+"_BD");
                }
                out.write(Consts.ls);
                int i = 0;
                for (String cell : cells) {
                    String cellToStore = cell.replaceAll("\"", "");
                    headerToIndex.put(cellToStore, i);
                    //System.out.println(cellToStore);
                    i++;
                }
            } else {
                //we process only the TCRB
                if (cells[headerToIndex.get("v_call")].equals("TRA")){
                    continue;
                }
                //in the case of CDR3+Vgene, DO NOT FORGET to convert to Immunoseq nomenclature
                String v = getVFamily(cells[headerToIndex.get("v_call")]); //we take the V family, EASIER to check
                String toCheck1 = cells[headerToIndex.get("junction_aa")];
                String cellIndex =cells[0].split("_")[1]; //cells[0].split("_")[0]; changed for the airr_cleaned file
                if (!mode.equals("cdr3")) {
                    toCheck1 += "-" + v;
                }
                //System.out.println(line);
                //go through each clonotypes
                for (String barcodeId : barcodeIdToTcr10X.keySet()) {
                        TCR tcr = barcodeIdToTcr10X.get(barcodeId);
                        String toCheck2 = tcr.getCdr3b();
                        if (!mode.equals("cdr3")) {
                            String v2 = getVFamily(tcr.getVb());
                            toCheck2 += "-" + v2;
                        }
                        //System.out.println("CHECK BD Rhaps "+toCheck1+ " with "+toCheck2);
                        if (toCheck1.equals(toCheck2)) {
                            System.out.println("--> HIT BD Cell_Index " + cellIndex + " " + toCheck1 + " with " + toCheck2);
                            tenxCells.add(tcr.getCellbarcode());
                            bdCells.add(cells[headerToIndex.get("cell_id")]);
                            //System.out.println("-->  SAME VDJ prot sequences !!");
                            int count = 1;
                            if (sequenceIdToHits.containsKey(barcodeId)) {
                                count += sequenceIdToHits.get(barcodeId);
                            }
                            sequenceIdToHits.put(barcodeId, count);
                            out.write(tcr.getSequence_id_B() + "\t" + tcr.getClone_id()+ "\t" + tcr.getClone_id_size()
                                    + "\t" + tcr.getDonor() + "\t" + tcr.getCell_type()
                                    + "\t" + tcr.getCellbarcode() + "\t" + tcr.getVb()
                                    + "\t" + tcr.getJb() + "\t" + tcr.getCdr3b()
                                    + "\t" + tcr.getVa()+ "\t" + tcr.getJa() + "\t" + tcr.getCdr3a() + "\t");
                            if (cellIndexToSarscov2.containsKey(tcr.getCellbarcode())) {
                                out.write(cellIndexToSarscov2.get(tcr.getCellbarcode()));
                            } else {
                                out.write("NA");
                            }
                            out.write("\t" + cellIndexToSample.get(cellIndex)
                                    + "\t" + cellIndexToCellType.get(cellIndex));
                            //write TCRA if there is
                            TCR tcrbd = barcodeIdToTcrBD.get(cells[headerToIndex.get("cell_id")]);
                            if (tcrbd.getCdr3a() !=null){
                                out.write("\t" +tcrbd.getVa()+"\t" +tcrbd.getJa()+"\t" +tcrbd.getCdr3a());
                            }
                            else{
                                out.write("\tNA\tNA\tNA");
                            }
                            out.write("\t" + line + Consts.ls);
                            if (tcr.getDonor().equals("BP_01")){
                                lineForSankeyGraphBP01.add("clone_"+ tcr.getClone_id()+" =["+ tcr.getClone_id_size()+","+
                                        cells[headerToIndex.get("clone_id")]+"]");
                            }
                            else if (tcr.getDonor().equals("BP_02")){
                                lineForSankeyGraphBP02.add("clone_"+ tcr.getClone_id()+" =["+ tcr.getClone_id_size()+","+
                                        cells[headerToIndex.get("clone_id")]+"]");
                            }


                            //we also check that TCRA matches too!
                            String tcra10x= tcr.getCdr3a()+"_"+ getVFamily(tcr.getVa());
                            //check if the tcr ahas a tcra
                            if (tcrbd.getCdr3a() != null) {
                                String tcrabd = tcrbd.getCdr3a() + "_" + getVFamily(tcrbd.getVa());
                                if (!tcrabd.equals(tcra10x)) {
                                    System.out.println("   --> NOOOOO TCRA hit " + tcra10x + " " + tcrabd);
                                    noTcrAhits++;
                                } else {
                                    System.out.println("        TCRA hit " + tcra10x + " " + tcrabd);
                                    tcrAhits++;
                                    int countA = 1;
                                    String bdBarcodeId= cells[headerToIndex.get("cell_id")];
                                    if (sequenceIdToHitsWithTRA.containsKey(bdBarcodeId)) {
                                        countA += sequenceIdToHitsWithTRA.get(bdBarcodeId);
                                    }
                                    sequenceIdToHitsWithTRA.put(bdBarcodeId, countA);
                                }
                            }
                            else{
                                missingTcrA++;
                                System.out.println("   --> NOOOOO TCRA hit because BD TCR has only TCRB!!! ");
                            }
                        }

                }
            }
            index++;
        }
        fileReader.close();
        out.close();
        System.out.println("number of cells from BD " + bdCells.size() + " are identical to " + tenxCells.size() + " 10X cells");
        System.out.println(" Number of cell with same TCRA too "+tcrAhits);
        System.out.println(" Number of cell without same TCRA "+noTcrAhits);
        System.out.println(" Number of cell with missing TCRA "+missingTcrA);
        System.out.println("number of hits found "+sequenceIdToHits.size());
        System.out.println("number of hits with TCRA found "+sequenceIdToHitsWithTRA.size());
        System.out.println("BP01");
        ArrayList<String> BP01 = new ArrayList<>(lineForSankeyGraphBP01);
        Collections.sort(BP01);
        for (String lin: BP01){
            System.out.println(lin);
        }
        ArrayList<String> BP02 = new ArrayList<>(lineForSankeyGraphBP02);
        Collections.sort(BP02);
        System.out.println("BP02");
        for (String lin: BP02){
            System.out.println(lin);
        }
    }


    private String getVFamily(String vGene) {
        String vFamily = vGene;
        if (vGene.contains("-")) {
            vFamily = vGene.split("-")[0];
        } else if (vGene.contains("*")) {
            vFamily = vGene.split("\\*")[0];
        }
        return vFamily;
    }


    public class TCR {
        private String sequence_id_B;
        private String sequence_id_A;
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
        private String cell_type;
        private String specificity = "NA"; //can be spike
        private int tra_hits = 0;
        private int trb_hits = 0;
        private String vdj_b_aa_seq;
        private String vdj_a_aa_seq;

        public TCR(String cellbarcode) {
            this.cellbarcode = cellbarcode;
        }

        public String getCellbarcode() {
            return cellbarcode;
        }

        public String getVdj_b_aa_seq() {
            return vdj_b_aa_seq;
        }

        public String getVdj_a_aa_seq() {
            return vdj_a_aa_seq;
        }

        public String getSequence_id_B() {
            return sequence_id_B;
        }

        public void setSequence_id_B(String sequence_id_B) {
            this.sequence_id_B = sequence_id_B;
        }

        public String getSequence_id_A() {
            return sequence_id_A;
        }

        public void setSequence_id_A(String sequence_id_A) {
            this.sequence_id_A = sequence_id_A;
        }

        public String getCell_type() {
            return cell_type;
        }

        public void setCell_type(String cell_type) {
            this.cell_type = cell_type;
        }

        public void setVdj_a_aa_seq(String vdj_a_aa_seq) {
            this.vdj_a_aa_seq = vdj_a_aa_seq;
        }

        public void setVdj_b_aa_seq(String vdj_aa_seq) {
            this.vdj_b_aa_seq = vdj_aa_seq;
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
