#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <fcntl.h>
#include <sys/stat.h>
#include <errno.h>
#include <string.h>
#include <time.h>

#include "headas_lock.h"

const char * HDtmpdir ()
{
	static char tmpdir[HD_TMPLEN];
	HDtmpdir_r(tmpdir);
	return tmpdir;
}


const char * HDtmpfile (const char * base, const char * ext)
{
	static char buffer[HD_TMPLEN];

	HDtmpfile_r(buffer, base, ext);

	return buffer;
}



int HDtmpdir_r (char * buffer)
{
	int code = 0;
	const char * tmpdir = getenv("HEADAS_TMPDIR");
	if (!tmpdir) {
		code = 1;
		tmpdir = ".";
	}
	strncpy(buffer, tmpdir, HD_TMPLEN);
	buffer[HD_TMPLEN-1]='\0';
	return code;
}



int HDtmpfile_r (char * buffer, const char * base, const char * ext)
{
	int code = 0;
	int fd;
	char tmpdir[HD_TMPLEN];
	char tmplock[HD_TMPLEN];
	struct flock lock = { 0 };

	buffer[0] = 0;

	if (!base)
		base = "";

	if (!ext)
		ext = ".tmp";

	HDtmpdir_r(tmpdir);

	sprintf(tmplock, "%s/%s", tmpdir, "lheatmp.lock");
	lock.l_type = F_WRLCK;

	fd = creat(tmplock, 0666);
	if (fd == -1) {
		printf("unable to open lock file '%s' [%s]\n",
			tmplock, strerror(errno));
		code = 1;
	}

	else if (fcntl(fd, F_SETLKW, &lock)) {
		printf("unable to lock temporary directory [%d, %s]\n",
				errno, strerror(errno));
		code = 2;
	}

	else {

		int probe = ((time(0) & 0xfff) << 20) | (getpid() << 4);

		while (1) {
			int test;
			struct stat st;
			sprintf(buffer, "%s/%s%x%s",
					tmpdir, base, probe, ext);
			test = stat(buffer, &st);
			if (test == -1) {
				if (errno == ENOENT)
					;
				else {
					code = 3;
					printf("error probing '%s' [%s]\n",
						buffer, strerror(errno));
				}
				break;
			}
			else
				++probe;
		}
	}

	close(fd);

	return code;
}
