package ch.chuv.bd;


/*
Copyright 2022 - Mathilde Foglierini Perez
This code is distributed open source under the terms of the GNU Free Documention License.
This class is used to convert the file '.....VDJ_Dominant_Contigs.csv' file to the usual 10X files used in changeo pipeline:
-filtered_contig.fasta
-filtered_contig_annotations.csv
It can be processed for TCR or BCR data, just change the variable isTCR. 
It can also select paired or non paired TCR (change the variable isPaired), because of the few paired TCRs we had with BD Rhapsody
 */

import ch.irb.ManageFastaFiles.FastaFileMaker;
import ch.irb.utils.Consts;

import java.io.*;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.LinkedHashMap;

public class BdRhapTo10Xconverter {

    private File folder = new File("/MYPATH/manuscript-tcell-fallet/processed_files/");
    private File outfolder = new File(folder.getPath()+Consts.fs+"TCR_10Xlike");
    private File bdFile = new File("/MYPATH/manuscript-tcell-fallet/BD_WTA_pipeline_out/BD_VDJ_Dominant_Contigs.csv");
    private HashMap<String,Integer> cellIndexToContigIndex = new HashMap<>();
    private ArrayList<String> selectedCellIndexes = new ArrayList<>();
    private HashMap<String, ArrayList<BdChain>> cellIndexToBdChains = new HashMap<>();

    //TODO HERE if BCR modify the code below
    public boolean isTCR=true;
    public boolean isPaired=false; //can be only TCRB to have more data (TCR mode only. for BCR we take only paired data)

    public static void main(String[] args) throws IOException {
        new BdRhapTo10Xconverter();
    }

    public BdRhapTo10Xconverter() throws IOException {

        if (!outfolder.exists()){
            outfolder.mkdir();
        }
        parseDbFile();
        write10XFiles();
    }

    private void parseDbFile() throws IOException {
        BufferedReader fileReader = new BufferedReader(new FileReader(bdFile));
        String line = "";
        int index = 0;
        int columnNumber = 0;
        int tcrNumb=0;
        int bcrNumb=0;
        int stopCodonSeq=0;
        int missingInfoSeq=0;
        HashMap<String, Integer> headerToIndex = new HashMap<>();
        while ((line = fileReader.readLine()) != null) {
            if (line.matches("#.*")) {
                continue;
            }

            String[] cells = line.split(",");
            if (index == 0) {
                //System.out.println("Get header "+line);
                int i = 0;
                for (String cell : cells) {
                    headerToIndex.put(cell, i);
                    i++;
                }
                columnNumber = i;
            } else {
                String cellIndex = cells[headerToIndex.get("Cell_Index")];
                String chainType = cells[headerToIndex.get("Chain_Type")];
                //we dont want TCR
                if (!isTCR & chainType.contains("TCR_")) {
                    tcrNumb++;
                    continue;
                } else if (isTCR & chainType.contains("BCR_")) {
                    bcrNumb++;
                    continue;
                }
                //we remove unproductive chain
                if (cells[headerToIndex.get("VDJ_Translation_Trimmed")].contains("*")) {
                    stopCodonSeq++;
                    continue;
                }
                //we remove sequence with missing info
                if (cells.length != columnNumber) {
                    missingInfoSeq++;
                    continue;
                }
                //System.out.println(line);
                BdChain bdChain = new BdChain(chainType, cellIndex,
                        cells[headerToIndex.get("C_gene_Dominant")], cells[headerToIndex.get("Read_Count")],
                        cells[headerToIndex.get("Molecule_Count")], cells[headerToIndex.get("V_gene_Dominant")],
                        cells[headerToIndex.get("D_gene_Dominant")], cells[headerToIndex.get("J_gene_Dominant")],
                        cells[headerToIndex.get("CDR3_Nucleotide_Dominant")], cells[headerToIndex.get("CDR3_Translation_Dominant")]);
                bdChain.setVdjNtSequence(cells[headerToIndex.get("Full_Contig_Nucleotide_Trimmed")]);
                ArrayList<BdChain> bdChains = new ArrayList<>();
                if (cellIndexToBdChains.containsKey(cellIndex)) {
                    bdChains = cellIndexToBdChains.get(cellIndex);
                }
                bdChains.add(bdChain);
                cellIndexToBdChains.put(cellIndex, bdChains);
            }
            index++;
        }
        fileReader.close();
        System.out.println("Number of removed TCR sequences: "+tcrNumb);
        System.out.println("Number of removed BCR sequences: "+bcrNumb);
        System.out.println("Number of removed codon stop sequences: "+stopCodonSeq);
        System.out.println("Number of removed missing info sequences: "+missingInfoSeq);
        System.out.println("Number of cells processed is " + cellIndexToBdChains.size());
        if (!isTCR) {
            for (String cellIndex : cellIndexToBdChains.keySet()) {
                ArrayList<BdChain> bdChains = cellIndexToBdChains.get(cellIndex);
                //System.out.println(cellIndex);
                boolean hasHeavy = false;
                boolean hasKappa = false;
                boolean hasLambda = false;
                double kappaUmis = 0;
                double lambdaUmis = 0;
                int kappaIndex = 0;
                int lambdaIndex = 0;
                int ind = 0;
                for (BdChain bdChain : bdChains) {
                    if (bdChain.getChainType().equals("IGH")) {
                        hasHeavy = true;
                    } else if (bdChain.getChainType().equals("IGK")) {
                        hasKappa = true;
                        kappaUmis = Double.parseDouble(bdChain.getUmis());
                        kappaIndex = ind;
                    } else if (bdChain.getChainType().equals("IGL")) {
                        hasLambda = true;
                        lambdaUmis = Double.parseDouble(bdChain.getUmis());
                        lambdaIndex = ind;
                    }
                    ind++;
                }
                //Here we check that we dont have an 'imbalanced' ratio between kappa and lambda chain, if yes we keep only 1
                if (hasKappa && hasLambda) {
                    boolean hasMoreKappa = true;
                    if (lambdaUmis > kappaUmis) {
                        hasMoreKappa = false;
                    }
                    if (hasMoreKappa) {//we remove the lambda chain
                        ArrayList<BdChain> chains = cellIndexToBdChains.get(cellIndex);
                        chains.remove(lambdaIndex);
                        cellIndexToBdChains.put(cellIndex, chains);
                    } else {//we remove the kappa chain
                        ArrayList<BdChain> chains = cellIndexToBdChains.get(cellIndex);
                        chains.remove(kappaIndex);
                        cellIndexToBdChains.put(cellIndex, chains);
                    }
                }

                // System.out.println(cellIndex+" has heavy "+hasHeavy+" has kappa "+hasKappa+" has lambda "+hasLambda);
                if (hasHeavy && (hasKappa || hasLambda)) {
                    selectedCellIndexes.add(cellIndex);
                }
            }
            System.out.println("After keeping only paired BCR data, Number of cells is " + selectedCellIndexes.size());
        }
        else{ //TCR case
            for (String cellIndex : cellIndexToBdChains.keySet()) {
                ArrayList<BdChain> bdChains = cellIndexToBdChains.get(cellIndex);
                //System.out.println(cellIndex);
                double betaUmis = 0;
                double alphaUmis = 0;
                double gammaUmis = 0;
                double deltaUmis = 0;
                boolean hasTCRB = false;
                boolean hasTCRA = false;
                boolean hasTCRG = false;
                boolean hasTCRD = false;
                int alphaIndex = 0;
                int betaIndex = 0;
                int gammaIndex = 0;
                int deltaIndex = 0;
                int ind = 0;
                for (BdChain bdChain : bdChains) {
                    //System.out.println("     "+bdChain.getChainType());
                    if (bdChain.getChainType().equals("TRB")) {
                        hasTCRB = true;
                        betaUmis = Double.parseDouble(bdChain.getUmis());
                        betaIndex = ind;
                    } else if (bdChain.getChainType().equals("TRA")) {
                        hasTCRA = true;
                        alphaUmis= Double.parseDouble(bdChain.getUmis());
                        alphaIndex = ind;
                    } else if (bdChain.getChainType().equals("TRG")) {
                        hasTCRG = true;
                        gammaUmis = Double.parseDouble(bdChain.getUmis());
                        gammaIndex = ind;
                    }else if (bdChain.getChainType().equals("TRD")) {
                        hasTCRD = true;
                        deltaIndex = ind;
                        deltaUmis  = Double.parseDouble(bdChain.getUmis());
                    }
                    ind++;
                }
                //Here we check that we dont have an 'imbalanced' ratio between beta and gamma chain, if yes we keep only 1
                if (hasTCRB && hasTCRD) {
                    boolean hasMoreBeta = true;
                    if (betaUmis < deltaUmis) {
                        hasMoreBeta = false;
                    }
                    if (hasMoreBeta) {//we remove the delta chain
                        ArrayList<BdChain> chains = cellIndexToBdChains.get(cellIndex);
                        chains.remove(deltaIndex);
                        cellIndexToBdChains.put(cellIndex, chains);
                        if (hasTCRG){ //and remove the gamma chain if exists
                            chains.remove(gammaIndex);
                            cellIndexToBdChains.put(cellIndex, chains);
                        }
                        hasTCRD=false;
                        hasTCRG=false;
                    } else {//we remove the beta chain
                        ArrayList<BdChain> chains = cellIndexToBdChains.get(cellIndex);
                        chains.remove(betaIndex);
                        cellIndexToBdChains.put(cellIndex, chains);
                        if (hasTCRA){ //and remove the alpha chain if exists
                            chains.remove(alphaIndex);
                            cellIndexToBdChains.put(cellIndex, chains);
                        }
                        hasTCRB=false;
                        hasTCRA=false;
                    }
                }

                //case where we still have 2 'light' chains
                if (hasTCRA && hasTCRG){
                    boolean hasMoreTCRA=true;
                    if (alphaUmis < gammaUmis){
                        hasMoreTCRA=false;
                    }
                    if(hasMoreTCRA){ //we remove the delta chain
                        ArrayList<BdChain> chains = cellIndexToBdChains.get(cellIndex);
                        chains.remove(gammaIndex);
                        cellIndexToBdChains.put(cellIndex, chains);
                        hasTCRG=false;
                    }else{ //we remove the alpha chain
                        ArrayList<BdChain> chains = cellIndexToBdChains.get(cellIndex);
                        chains.remove(alphaIndex);
                        cellIndexToBdChains.put(cellIndex, chains);
                        hasTCRA=false;
                    }
                }

                // System.out.println(cellIndex+" has heavy "+hasHeavy+" has kappa "+hasKappa+" has lambda "+hasLambda);
                if (isPaired) {
                    if ((hasTCRB && hasTCRA) || (hasTCRG && hasTCRD)) {
                        selectedCellIndexes.add(cellIndex);
                    }
                } else if (hasTCRB) { //for Fallet project we take the TCRB to have more data
                    selectedCellIndexes.add(cellIndex);
                }
            }
            System.out.println("After keeping only paired ("+isPaired+") TCR data, Number of cells is " + selectedCellIndexes.size());
        }

    }

    private void write10XFiles() throws IOException {
        LinkedHashMap<String,String> idToSeq = new LinkedHashMap<>();
        File outFile = new File(outfolder.getPath()+Consts.fs+"filtered_contig_annotations.csv");
        BufferedWriter out = new BufferedWriter(new FileWriter(outFile));
        out.write("barcode,contig_id,locus,c_gene,reads,umis,v_gene,d_gene,j_gene,cdr3_nt,cdr3"+Consts.ls);
        for (String cellIndex: selectedCellIndexes){
            ArrayList<BdChain> bdChains = cellIndexToBdChains.get(cellIndex);
            for(BdChain bdChain: bdChains){
                out.write(bdChain.getCellIndex()+","+bdChain.getContigId()+","+bdChain.getChainType()+","
                +bdChain.getcGene()+","+bdChain.getReads()+","+bdChain.getUmis()+","+bdChain.getvGene()+","
                + bdChain.getdGene()+","+bdChain.getjGene() +","+bdChain.getCdr3Nt()+","+bdChain.getCdr3aa()
                + Consts.ls);
                idToSeq.put(bdChain.getContigId(), bdChain.getVdjNtSequence());
            }
        }
        out.close();
        //write fasta file
        new FastaFileMaker(outfolder.getPath()+Consts.fs+"filtered_contig.fasta",idToSeq);
        System.out.println("Number of fasta id written "+idToSeq.size());
    }

    private class BdChain{
        private String chainType;
        private String cellIndex;
        private String cGene;
        private String reads;
        private String umis;
        private String vGene;
        private String dGene;
        private String jGene;
        private String cdr3Nt;
        private String cdr3aa;
        private String contigId;
        private String vdjNtSequence;

        public BdChain(String chainType, String cellIndex, String cGene, String reads, String umis, String vGene, String dGene, String jGene, String cdr3Nt, String cdr3aa) {
            this.chainType = chainType;
            this.cellIndex = cellIndex;
            this.cGene = cGene;
            this.reads = reads;
            this.umis = umis;
            this.vGene = vGene;
            this.dGene = dGene;
            this.jGene = jGene;
            this.cdr3Nt = cdr3Nt;
            this.cdr3aa = cdr3aa;
            setContigId();
        }

        public String getContigId() {
            return contigId;
        }

        public void setContigId() {
            int contigIndex=1;
            if (cellIndexToContigIndex.containsKey(cellIndex)){
                contigIndex+= cellIndexToContigIndex.get(cellIndex);
            }
            this.contigId = cellIndex+"_contig_"+contigIndex;
            cellIndexToContigIndex.put(cellIndex,contigIndex);
        }

        public String getVdjNtSequence() {
            return vdjNtSequence;
        }

        public void setVdjNtSequence(String vdjNtSeq) {
            this.vdjNtSequence = vdjNtSeq;
        }

        public String getChainType() {
            String locus="IGH";
            if (chainType.equals("BCR_Kappa")){
                locus="IGK";
            } else if (chainType.equals("BCR_Lambda")){
                locus="IGL";
            } else if (chainType.equals("TCR_Alpha")) {
                locus="TRA";
            } else if (chainType.equals("TCR_Beta")) {
                locus="TRB";
            } else if (chainType.equals("TCR_Gamma")) {
                locus="TRG";
            } else if (chainType.equals("TCR_Delta")) {
                locus="TRD";
            }
            return locus;
        }

        public String getCellIndex() {
            return cellIndex;
        }

        public String getcGene() {
            return cGene;
        }

        public String getReads() {
            return reads;
        }

        public String getUmis() {
            return umis;
        }

        public String getvGene() {
            return vGene;
        }

        public String getdGene() {
            return dGene;
        }

        public String getjGene() {
            return jGene;
        }

        public String getCdr3Nt() {
            return cdr3Nt;
        }

        public String getCdr3aa() {
            return cdr3aa;
        }
    }

}
