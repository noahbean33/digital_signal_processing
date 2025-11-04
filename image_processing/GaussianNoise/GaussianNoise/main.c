#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define BMP_HEADER_SIZE             54
#define BMP_COLOR_TABLE_SIZE        1024
#define CUSTOM_IMG_SIZE             512*512


void imageWriter(const char *imgName,
                 unsigned char *header,
                 unsigned char *colorTable,
                 unsigned char *buf,
                 int bitDepth);

void  imageReader(const char *imgName,
                  int *_height,
                  int *_width,
                  int *_bitDepth,
                  unsigned char *_header,
                  unsigned char *_colorTable,
                  unsigned char *_buf
                  );
void gaussian(unsigned char *_inputImgData, int imgCols, int imgRows, float var, float mean);

int main()
{
    int imgWidth,imgHeight,imgBitDepth;
    unsigned char imgHeader[BMP_HEADER_SIZE];
    unsigned char imgColorTable[BMP_COLOR_TABLE_SIZE];
    unsigned char imgBuffer[CUSTOM_IMG_SIZE];
    unsigned char imgBuffer2[CUSTOM_IMG_SIZE];

    const char imgName[] = "lena512.bmp";
    const char newImgName[] = "lena_gauss.bmp";

    imageReader(imgName,&imgWidth,&imgHeight,&imgBitDepth,&imgHeader[0],&imgColorTable[0],&imgBuffer[0]);
    gaussian(&imgBuffer[0],imgWidth,imgHeight,90,80);

    imageWriter(newImgName,imgHeader,imgColorTable,imgBuffer,imgBitDepth);
    return 0;
}

void imageWriter(const char *imgName,
                 unsigned char *header,
                 unsigned char *colorTable,
                 unsigned char *buf,
                 int bitDepth)
   {
     FILE *fo = fopen(imgName,"wb");
     fwrite(header,sizeof(unsigned char),54,fo);
     if(bitDepth <=8)
     {
         fwrite(colorTable,sizeof(unsigned char),1024,fo);
     }
     fwrite(buf,sizeof(unsigned char),CUSTOM_IMG_SIZE, fo);
     fclose(fo);

   }


   void  imageReader(const char *imgName,
                  int *_height,
                  int *_width,
                  int *_bitDepth,
                  unsigned char *_header,
                  unsigned char *_colorTable,
                  unsigned char *_buf
                  )
{
    int i;
    FILE *streamIn;
    streamIn = fopen(imgName,"rb");

    if(streamIn ==(FILE *)0)
    {

        printf("Unable to read image \n");
    }

    for(i =0;i<54;i++)
    {
        _header[i] = getc(streamIn);
    }

    *_width = *(int *)&_header[18];
    *_height = *(int *)&_header[22];
    *_bitDepth = *(int *)&_header[28];

    if(*_bitDepth <=8)
    {
        fread(_colorTable,sizeof(unsigned char),1024,streamIn);
    }

    fread(_buf,sizeof(unsigned char),CUSTOM_IMG_SIZE,streamIn);

    fclose(streamIn);
}

void gaussian(unsigned char *_inputImgData, int imgCols, int imgRows, float var, float mean)
{
    int x,y;
    float noise, theta;
    for(y =0;y<imgRows;y++)
        for(x=0;x<imgCols;x++)
    {
        noise = sqrt(-2*var*log(1.0-(float)rand()/32767.1));
        theta = rand()*1.9175345e-4 - 3.14159265;
        noise = noise * cos(theta);
        noise = noise + mean;
        if(noise >255)noise = 255;
        if(noise <0) noise = 0;
       *(_inputImgData+x+(long)y*imgCols) = (unsigned char)(noise +0.5);
    }
}
