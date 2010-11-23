#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <fcntl.h>
#include <sys/mman.h>

#define FILEPATH "/tmp/mmapped.bin"
#define NUMINTS  (1000)
#define FILESIZE (NUMINTS * sizeof(int))

int main(int argc, char *argv[])
{
  int i;
  int fd;
  int *map;  /* mmapped array of int's */

  (void)argc;
  (void)argv;

  fd = open(FILEPATH, O_RDONLY);
  if (fd == -1) {
    perror("Error opening file for reading");
    exit(EXIT_FAILURE);
  }

  map = (int*)mmap(0, FILESIZE, PROT_READ, MAP_SHARED, fd, 0);
  if (map == MAP_FAILED) {
    close(fd);
    perror("Error mmapping the file");
    exit(EXIT_FAILURE);
  }
    
  /* Read the file int-by-int from the mmap
   */
  for (i = 0; i <NUMINTS; ++i) {
    printf("%d: %d\n", i, map[i]);
  }

  if (munmap(map, FILESIZE) == -1) {
    perror("Error un-mmapping the file");
  }
  close(fd);
  return 0;
}
