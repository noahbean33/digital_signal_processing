#include <stdio.h>
#include <stdlib.h>

int main()
{
   FILE * fp = fopen("girlface.bmp","rb");
  unsigned char  *imageData;
  unsigned char  *newImageData;
  unsigned char imageHeader[54];
  unsigned char colorTable[1024];

  fread(imageHeader,sizeof(unsigned char),54,fp);

  int width = *(int*)&imageHeader[18];
  int height = *(int*)&imageHeader[22];
  int bitDepth = *(int*)&imageHeader[28];

  int imgDataSize = width * height;

  imageData = (unsigned char *)malloc(imgDataSize *sizeof(unsigned char));
  newImageData = (unsigned char *)malloc(imgDataSize *sizeof(unsigned char));

  if(bitDepth <=8)
  {
      fread(colorTable,sizeof(unsigned char),1024,fp);
  }
     fread(imageData,sizeof(unsigned char),imgDataSize,fp);

     //Computing image negative
     for(int i =0; i<height;i++)
     {
         for(int j =0;j<width;j++)
         {
             newImageData[i*width+j] =  255 - imageData[i*width+j];
         }

     }

     FILE *fo = fopen("girlface_negative.bmp","wb");
     fwrite(imageHeader, sizeof(unsigned char),54,fo);

     if(bitDepth <=8)
     {
         fwrite(colorTable,sizeof(unsigned char),1024,fo);
     }

     fwrite(newImageData,sizeof(unsigned char),imgDataSize,fo);

     printf("Success !\n");

     fclose(fo);
     fclose(fp);

    return 0;
}
