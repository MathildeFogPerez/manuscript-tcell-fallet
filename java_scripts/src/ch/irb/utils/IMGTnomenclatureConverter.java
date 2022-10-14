package ch.irb.utils;


/**
 * Created by Mathilde on 03.04.2017.
 * This class is used to get the IMGT (ex: Homsap TRBV6-5*01) or ImmunoSeq nomenclature (ex: TCRBV06-05*01)
 */
public class IMGTnomenclatureConverter {

    private String imgtGene;
    private String immunoSeqGene;

    public static void main(String[] args) {
        IMGTnomenclatureConverter converter = new IMGTnomenclatureConverter("Homsap TRBV7-8*01 F", true);
        System.out.println(converter.getImgtGene()+" to ImmunoSeq "+converter.getImmunoSeqGene());
    }

    public IMGTnomenclatureConverter(String gene, boolean isIMGT) {
        String geneToProcess ="V";
        if (gene.contains("J")){
            geneToProcess = "J";
        }
        //System.out.println(gene);
        gene=gene.trim();
        if (isIMGT) {
            imgtGene = gene;
            String modifSeqGene = gene.replace("Homsap ", "").replace(" (F)","").replace("[F]","").replace(" F","").replace("TR", "TCR").trim();
            String firstNumber = modifSeqGene.split("[VJ]")[1].split("\\*")[0]; //there is always the allele with IMGT
            String modifiedFirstNumber = firstNumber;
            String modifiedSecNumber = "";
            String allele = "";
            if (gene.contains("-")){
                modifiedFirstNumber = firstNumber.split("-")[0];
                String secondNumber = modifSeqGene.split("[VJ]")[1].split("-")[1];
                //System.out.println("Second number is "+secondNumber);
                if (secondNumber.contains("*")) { //there is the allele
                    secondNumber = secondNumber.split("\\*")[0];
                }
                modifiedSecNumber = "-"+secondNumber;
                if (secondNumber.matches("\\d{1}")) {
                    modifiedSecNumber = "-0" + secondNumber;
                }
                //System.out.println("Modified Second number is "+modifiedSecNumber);
            }else if (!gene.contains("TRBV26")){
                modifiedSecNumber="-01";
            }
            if (modifSeqGene.contains("*")) { //there is the allele
                allele = "*" + modifSeqGene.split("\\*")[1];
            }
            if (modifiedFirstNumber.matches("\\d{1}")) {
                modifiedFirstNumber = "0" + modifiedFirstNumber;
            }
            String immunoSeqvGene = modifSeqGene.split("[VJ]")[0] +geneToProcess+ modifiedFirstNumber + modifiedSecNumber;
            String modifImmunoSeqvGene= checkIfVGeneExist(immunoSeqvGene);
            immunoSeqGene = modifImmunoSeqvGene + allele;
            //System.out.println("IMGT: " + imgtGene + " to ImmunoSeq " + immunoSeqGene);
        } else {
            //System.out.println(gene);
            immunoSeqGene = gene;
            if (!gene.contains("TCR")){//case where we have the nuc sequence of the rearrangement
                imgtGene="notKnown";
                return;
            }
            String modifSeqGene = gene.replace("TCR", "TR");
            String firstNumber = modifSeqGene.split("[VJ]")[1].split("-")[0];
            String modifiedFirstNumber = firstNumber;
            if (firstNumber.matches("0\\d{1}")) {
                modifiedFirstNumber = firstNumber.replaceFirst("0", "");
            }
            String modifiedSecNumber = "";
            String allele = "";
            if (gene.contains("-")) {
                String secondNumber = modifSeqGene.split("[VJ]")[1].split("-")[1];
                allele = "";
                if (secondNumber.contains("*")) { //there is the allele
                    secondNumber = secondNumber.split("\\*")[0];
                    allele = "*" + modifSeqGene.split("\\*")[1];
                }
                modifiedSecNumber = "-" +secondNumber;
                if (secondNumber.matches("0\\d{1}")) {
                    modifiedSecNumber ="-" + secondNumber.replaceFirst("0", "");
                }
            }
            imgtGene = modifSeqGene.split("[VJ]")[0]+geneToProcess + modifiedFirstNumber + modifiedSecNumber + allele;
            //System.out.println("ImmunoSeq: " + immunoSeqGene + " to IMGT " + imgtGene);
        }
    }

    public String getImgtGene() {
        return imgtGene;
    }

    public void setImgtGene(String imgtGene) {
        this.imgtGene = imgtGene;
    }

    public String getImmunoSeqGene() {
        return immunoSeqGene;
    }

    public String getImmunoSeqFamily(){
        if (immunoSeqGene.contains("-")){
            return immunoSeqGene.split("-")[0];
        }
        else {
            return immunoSeqGene;
        }
    }

    public void setImmunoSeqGene(String immunoSeqGene) {
        this.immunoSeqGene = immunoSeqGene;
    }

    private String checkIfVGeneExist(String immunoSeqVgene){
        String vGene = immunoSeqVgene;
        //check if exist in the list
        //has to be updated below because it is old and does to include v gene!
        /*if (!Consts.immunoSeqVGenes.contains(vGene)){
            //if not we return the family
           // System.out.println("We dont have the vGENE "+vGene+" we take the vFAMILY "+immunoSeqVgene.split("-")[0]);
            vGene = immunoSeqVgene.split("-")[0];
        }*/
        return vGene;
    }
}
