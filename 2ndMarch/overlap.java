package overlap;

import java.awt.Color;
import java.awt.image.BufferedImage;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;

import javax.imageio.ImageIO;

class vector {
    int minx = 0;
    int miny = 0;
    int maxx = 0;
    int maxy = 0;
}

public class overlap {

    //vector data
    static vector[] vecotrs = new vector[100];
    static vector[] vecotrs2 = new vector[100];

    static int[] area1 = new int[100];
    static int[] area2 = new int[100];
    static int[] overlap = new int[100];
    static int im1OB = 0;
    static int im2OB = 0;

    public static Color[][] loadPixelsFromImage(File file, File file2) throws IOException {

        //read the image1 and image2 file
        BufferedImage image = ImageIO.read(file);
        BufferedImage image2 = ImageIO.read(file2);

        //this is for the out put image (result)
        BufferedImage im = new BufferedImage(image.getWidth(), image.getHeight(), BufferedImage.TYPE_INT_RGB);

        //make the 2d array for image1  and image2 
        Color[][] colors = new Color[image.getWidth()][image.getHeight()];
        Color[][] colors2 = new Color[image2.getWidth()][image2.getHeight()];

        //make the 2d array for changed detected pixels
        Color[][] change = new Color[image.getWidth()][image.getHeight()];

        //init for debug
        for (int x = 0; x < image2.getWidth(); x++) {
            for (int y = 0; y < image2.getHeight(); y++) {
                change[x][y] = new Color(0, 0, 0);
                colors[x][y] = new Color(0, 0, 0);
                colors2[x][y] = new Color(0, 0, 0);
            }
        }

        /*making the loop for the 3 objects for image1*/
        for (int i = 0; i < im1OB; i++) {
            //get image1 pixel data 
            for (int x = vecotrs[i].minx; x < vecotrs[i].maxx; x++) {
                for (int y = vecotrs[i].miny; y < vecotrs[i].maxy; y++) {
                    colors[x][y] = new Color(image.getRGB(x, y));
                    //object1
                    if (colors[x][y].getBlue() > 100) {
                        area1[i]++;
                    }

                }
            }
        }
        for (int i = 0; i < im2OB; i++) {
            //get image1 pixel data 
            for (int x = vecotrs2[i].minx; x < vecotrs2[i].maxx; x++) {
                for (int y = vecotrs2[i].miny; y < vecotrs2[i].maxy; y++) {
                    colors2[x][y] = new Color(image2.getRGB(x, y));
                    if (colors2[x][y].getBlue() > 100) {
                        area2[i]++;
                    }
                }
            }
        }

        for (int i = 0; i < im1OB; i++) {
            for (int x = vecotrs[i].minx; x < vecotrs[i].maxx; x++) {
                for (int y = vecotrs[i].miny; y < vecotrs[i].maxy; y++) {
                    // if(x>220 && x<315 && y>60 &&y<200){
                    int red = Math.abs(colors[x][y].getRed() + colors2[x][y].getRed());
                    int blue = Math.abs(colors[x][y].getBlue() + colors2[x][y].getBlue());
                    int green = Math.abs(colors[x][y].getGreen() + colors2[x][y].getGreen());

                    //out boundary
                    if (red < 30 && blue < 30 && green < 30)
                        change[x][y] = new Color(0, 0, 0);
                    else if (red < 450 && blue < 450 && green < 450)//change 
                        change[x][y] = new Color(0, 0, 255);
                    else if (red >= 450 && blue >= 450 && green >= 450) {//overlap
                        change[x][y] = new Color(255, 255, 255);
                        overlap[i]++;
                    }
                }
            }

        }

        for (int i = 0; i < image.getWidth(); i++)
            for (int j = 0; j < image.getHeight(); j++)
                im.setRGB(i, j, change[i][j].getRGB());

        ImageIO.write(im, "PNG", new File("Image1.jpg"));

        return colors;
    }

    public static int readfile(String text, int file) {
        //read file data
        FileReader fr = null;
        BufferedReader br = null;
        int jj = 0;
        String ob = "ob" + 1;
        try {
            fr = new FileReader(text);
            br = new BufferedReader(fr);

            String line;
            int MINX = 1000, MAXX = 0, MINY = 1000, MAXY = 0;
            while ((line = br.readLine()) != null) {
            
                if (line.equals(ob)) {    //object list/ object is separated by "ob+number"
                    jj++;
                    ob = "ob" + (jj + 1);
                    MINX = 1000;    MAXX = 0;
                    MINY = 1000;    MAXY = 0;
                } else {
                    //token 0 = x, token1 =y
                    String[] tokens = line.split(" ");
                    int x = Integer.parseInt(tokens[0]);
                    int y = Integer.parseInt(tokens[1]);
                    //get image objects 
                    if (file == 1) {//image 1
                        if (MINX > x) {
                            vecotrs[jj].minx = x;
                            MINX = x;
                        }
                        if (MAXX < x) {
                            vecotrs[jj].maxx = x;
                            MAXX = x;
                        }
                        if (MINY > y) {
                            vecotrs[jj].miny = y;
                            MINY = y;
                        }
                        if (MAXY < y) {
                            vecotrs[jj].maxy = y;
                            MAXY = y;
                        }
                    } else {// image 2
                        if (MINX > x) {
                            vecotrs2[jj].minx = x;
                            MINX = x;
                        }
                        if (MAXX < x) {
                            vecotrs2[jj].maxx = x;
                            MAXX = x;
                        }
                        if (MINY > y) {
                            vecotrs2[jj].miny = y;
                            MINY = y;
                        }
                        if (MAXY < y) {
                            vecotrs2[jj].maxy = y;
                            MAXY = y;
                        }
                    }
                }
            }
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        } catch (IOException e) {
            e.printStackTrace();
        } finally {
            try {
                br.close();
                fr.close();
            } catch (IOException e) {
                e.printStackTrace();
            }
        }
        return jj;

    }

    public static void main(String[] args) throws IOException {

        //init vec
        for (int i = 0; i < 100; i++) {
            vecotrs[i] = new vector();
            vecotrs2[i] = new vector();
        }

        im1OB = readfile("t.txt", 1);
        im2OB = readfile("t2.txt", 2);

        //chnage detection
        Color[][] colors = loadPixelsFromImage(new File("t5.png"), new File("t6.png"));

        
        // display result
        for (int i = 0; i < im1OB; i++) {
            double over = ((double) overlap[i] / (double) area1[i]) * 100;

            System.out.println("The area1 area = " + area1[i]);
            System.out.println("The area2 area = " + area2[i]);
            System.out.println("The overlap area = " + overlap[i]);

            System.out.println("The percentage of overlap area:\n " + overlap[i] + "/" + area1[i] + " = " + (int) over + "%\n\n\n");
        }

    }
}
