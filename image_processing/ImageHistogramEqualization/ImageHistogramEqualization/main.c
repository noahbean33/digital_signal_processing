#include <stdio.h>
#include <stdlib.h>
#define BMP_HEADER_SIZE         54
#define BMP_COLOR_TABLE_SIZE    1024
#define CUSTOM_IMG_SIZE         512*512

float IMG_HIST[255];

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
                 int bitDepth) ;

void ImgHistogram(unsigned char * _imgData, int imgRows, int imgCols, float hist[]);
void ImgHistogramEqualization(unsigned char *_inputImgData, unsigned char *_outputImgData,int imgRows,int imgCols);

int main()
{
    int imgWidth, imgHeight,imgBitDepth;
    unsigned char imgHeader[BMP_HEADER_SIZE];
    unsigned char imgColorTable[BMP_COLOR_TABLE_SIZE];
    unsigned char imgBuffer[CUSTOM_IMG_SIZE];
    unsigned char imgBuffer2[CUSTOM_IMG_SIZE];

    const char imgName[] = "lena512.bmp";
    const char newImgName[] ="lena_eqz.bmp";

    imageReader(imgName,&imgHeight,&imgWidth,&imgBitDepth,&imgHeader[0],&imgColorTable[0],&imgBuffer[0]);
    ImgHistogramEqualization(&imgBuffer[0],&imgBuffer2[0],imgHeight,imgWidth);
    imageWriter(newImgName,imgHeader,imgColorTable,imgBuffer2,imgBitDepth);

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


void ImgHistogram(unsigned char * _imgData, int imgRows, int imgCols, float hist[])
{
  FILE *fptr;
  fptr = fopen("image_hist.txt","w");
  int x,y,i,j;
  long int ihist[255];
  long int sum;

  for(i=0;i<=255;i++)
  {
      ihist[i] = 0;
  }
  sum =0;
  for(y =0;y<imgRows;y++)
  {
      for(x=0;x<imgCols;x++)
      {
          j = *(_imgData+x+y*imgCols);
          ihist[j] = ihist[j] +1;
          sum = sum +1;

      }
  }

    for( i =0;i<255;i++)
    {
        hist[i] = (float)ihist[i]/(float)sum;

    }


     for(int i=0;i<255;i++)
     {
         fprintf(fptr,"\n%f",hist[i]);
     }

    fclose(fptr);

}

void ImgHistogramEqualization(unsigned char *_inputImgData, unsigned char *_outputImgData,int imgRows,int imgCols)
{
    int x,y,i,j;
    int histeq[256];
    float hist[256];
    float sum;

    ImgHistogram(&_inputImgData[0], imgRows,imgCols,&hist[0]);

    for(i=0;i<255;i++)
    {
        sum =0.0;
        for(j=0;j<=i;j++)
        {
            sum = sum+hist[j];
        }
        histeq[i] = (int)(255*sum+0.5);
    }
    for(y =0;y<imgRows;y++)
    {
        for(x=0;x<imgCols;x++)
        {
            *(_outputImgData+x+y*imgCols) =  histeq[*(_inputImgData+x+y*imgCols)];
        }
    }
}





















