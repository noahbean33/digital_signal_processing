#include <stdio.h>
#include <stdlib.h>

#define BMP_HEADER_SIZE         54
#define BMP_COLOR_TABLE_SIZE    1024
#define CUSTOM_IMG_SIZE         512*512

void median_filter(unsigned char *_inputImgData, unsigned char *_outputImgData, int imgCols,int imgRows);
    void  imageReader(const char *imgName,
                  int *_height,
                  int *_width,
                  int *_bitDepth,
                  unsigned char *_header,
                  unsigned char *_colorTable,
                  unsigned char *_buf
                  );

    void imageWriter(const char *imgName,
                 unsigned char *header,
                 unsigned char *colorTable,
                 unsigned char *buf,
                 int bitDepth);
int main()
{
    int imgWidth, imgHeight,imgBitDepth;
    unsigned char imgHeader[BMP_HEADER_SIZE];
    unsigned char imgColorTable[BMP_COLOR_TABLE_SIZE];
    unsigned char imgBuffer[CUSTOM_IMG_SIZE];
    unsigned char imgBuffer2[CUSTOM_IMG_SIZE];

    const char imgName[] ="lena512.bmp";
    const char newImgName[] = "lena_med_filt.bmp";

    imageReader(imgName,&imgHeight,&imgWidth,&imgBitDepth,&imgHeader[0],&imgColorTable[0],&imgBuffer[0]);
    median_filter(&imgBuffer[0],&imgBuffer2[0],imgWidth,imgHeight);

    imageWriter(newImgName,imgHeader,imgColorTable,&imgBuffer2[0],imgBitDepth);
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

void median_filter(unsigned char *_inputImgData,unsigned char *_outputImgData,int imgCols,int imgRows)
{
   int x,y,i,j,z;
   int n, ar[121],a;
   n =7;

   for(y = n/2; y<imgRows-n/2;y++)
    for(x = n/2;x<imgCols-n/2;x++)
   {
       z=0;
       for(j=-n/2;j<=n/2;j++)
       for(i=-n/2;i<=n/2;i++){
       ar[z] = *(_inputImgData+x+i+(long)(y+j)*imgCols);
       z++;

   }
   for(j =1;j<=n*n-1;j++)
   {

       a = ar[j];
       i= j-1;
       while(i>=0&&ar[i] > a)
       {

           ar[i+1] = ar[i];
           i =i-1;
       }
       ar[i+1] =a;
   }
   *(_outputImgData +x+(long)y*imgCols )= ar[n*n/2];
   }
}
