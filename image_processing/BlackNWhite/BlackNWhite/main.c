#include <stdio.h>
#include <stdlib.h>

#define WHITE   255
#define BLACK   0
#define THRESHOLD   50


int main()
{

   FILE *fIn = fopen("images/lighthouse.bmp","rb");
   FILE *fOut = fopen("images/lighthouse_bw2.bmp","wb");

   unsigned  char imgHeader[54];
   unsigned  char colorTable[1024];

   if(fIn == NULL)
   {
       printf("Unable to open image \n");
   }

   for(int i=0;i<54;i++)
   {
       imgHeader[i] = getc(fIn);
   }

   fwrite(imgHeader,sizeof(unsigned char),54,fOut);

   int height = *(int*)&imgHeader[22];
   int width  = *(int*)&imgHeader[18];
   int bitDepth = *(int *)&imgHeader[28];

   int imgSize = height * width;

   if(bitDepth<=8)
   {
       fread(colorTable, sizeof(unsigned char),1024,fIn);
       fwrite(colorTable,sizeof(unsigned char),1024,fOut);
   }

   unsigned char buffer[imgSize];
   fread(buffer,sizeof(unsigned char),imgSize,fIn);

    //Black and White converter
    for(int i =0;i<imgSize;i++)
    {
        buffer[i] = (buffer[i]>THRESHOLD)? WHITE : BLACK;
    }

    fwrite(buffer,sizeof(unsigned char), imgSize,fOut);
    fclose(fIn);
    fclose(fOut);
    return 0;
}
