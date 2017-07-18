import java.io.FileWriter;
import java.io.IOException;

public class OutPut {
    public static void writeOut(String outfile, String toWrite) throws IOException{
    //String fileName = "C:/Users/venkav3/Dropbox/Research/new/af/VFRP15/Algorithm/Data/Xtrial.txt";
    //String toWrite = "message to write " + Double.toString(test2) + " " + Integer.toString(test1) + "\n";
    
    	
    FileWriter out1 = new FileWriter(outfile, true);        // will overwrite existing file, if true, it will append
    out1.write(toWrite);
    out1.close();
    }
    public static void initFile(String outfile) throws IOException{
        //String fileName = "C:/Users/venkav3/Dropbox/Research/new/af/VFRP15/Algorithm/Data/Xtrial.txt";
        //String toWrite = "message to write " + Double.toString(test2) + " " + Integer.toString(test1) + "\n";
        
        	
        FileWriter out1 = new FileWriter(outfile, false);        // will overwrite existing file, if true, it will append
        out1.write("\r\n");
        out1.close();
        }
}

