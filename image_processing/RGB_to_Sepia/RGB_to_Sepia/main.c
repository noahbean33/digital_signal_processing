#include <stdio.h>
#include <stdlib.h>

#define MAX_PIXEL   255

int main()
{
    FILE  *fIn =fopen("images/lena_color.bmp","rb");
    FILE  *fOut =fopen("images/lena_sepia.bmp","wb");

    unsigned char imgHeader[54];
    unsigned char colorTable[1024];

    if(fIn ==NULL)
    {
        printf("Unable to open file \n");
    }
    for(int i =0;i<54;i++)
    {
       imgHeader[i] = getc(fIn);
    }
    fwrite(imgHeader,sizeof(unsigned char),54,fOut);

    int height = *(int*)&imgHeader[22];
    int width  = *(int*)&imgHeader[18];
    int bitDepth = *(int*)&imgHeader[28];

    if(bitDepth <=8)
    {
        fread(colorTable,sizeof(unsigned char),1024,fIn);
        fwrite(colorTable,sizeof(unsigned char),1024,fOut);
    }

    int imgSize = height * width;
    int r,g,b;
    unsigned char buffer[imgSize][3];

    for(int i =0;i<imgSize;i++)
    {
        r=g=b=0;
        buffer[i][0] = getc(fIn);      //red
        buffer[i][1] = getc(fIn);      //green
        buffer[i][2] = getc(fIn);      //blue

        r = (buffer[i][0]*0.393) + (buffer[i][1]*0.769) + (buffer[i][2]*0.189);
        g = (buffer[i][0]*0.349) + (buffer[i][1]*0.686) + (buffer[i][2]*0.168);
        b = (buffer[i][0]*0.272) + (buffer[i][1]*0.534) + (buffer[i][2]*0.131);

        if(r>MAX_PIXEL){ r =MAX_PIXEL;}
        if(g>MAX_PIXEL){ g = MAX_PIXEL;}
        if(b>MAX_PIXEL){ b = MAX_PIXEL;}

        putc(b,fOut);
        putc(g,fOut);
        putc(r,fOut);

    }

    fclose(fOut);
    fclose(fIn);

    return 0;
}
