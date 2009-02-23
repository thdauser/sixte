/* print_errno_info.c - A function that sends an appropriate error message to stderr
 *
 * Copyright (C) 2005 - 2006
 *                    Joern Wilms
 *                    Robert Dowse
 *                    Robert Swain
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
 *
 */

#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <sysexits.h>


void print_errno_info() { //int errnumber) {
  fprintf(stderr,"Function returned errno=%i\n",errno);

  switch(errno) {
  // complete list valid for Linux, kernel 2.4

  // 124 EMEDIUMTYPE      Wrong medium type
#ifdef EMEDIUMTYPE
 case EMEDIUMTYPE:
   fprintf(stderr,"124 EMEDIUMTYPE: Wrong medium type\n");
   break;
#endif
   // 123 ENOMEDIUM        No medium found
#ifdef ENOMEDIUM
 case ENOMEDIUM:
   fprintf(stderr,"123 ENOMEDIUM: No medium found\n");
#endif
   // 122 EDQUOT           Disk quota exceeded
#ifdef EDQUOT
 case EDQUOT:
   fprintf(stderr,"122 EDQUOT: Disk quota exceeded\n");
   break;
#endif
   // 121 EREMOTEIO        Remote I/O error
#ifdef EREMOTEIO
 case EREMOTEIO:
   fprintf(stderr,"121 EREMOTEIO: Remote I/O error\n");
   break;
#endif
   // 120 EISNAM           Is a named type file
#ifdef EISNAM
 case EISNAM:
   fprintf(stderr,"120 EISNAM: Is a named type file\n");
   break;
#endif
   // 119 ENAVAIL          No XENIX semaphores available
#ifdef ENAVAIL
 case ENAVAIL:
   fprintf(stderr,"119 ENAVAIL: No XENIX semaphores available\n");
   break;
#endif
   // 118 ENOTNAM          Not a XENIX named type file
#ifdef ENOTNAM
 case ENOTNAM:
   fprintf(stderr,"118 ENOTNAM: Not a XENIX named type file\n");
   break;
#endif
   // 117 EUCLEAN          Structure needs cleaning
#ifdef EUCLEAN
 case EUCLEAN:
   fprintf(stderr,"117 EUCLEAN: Structure needs cleaning\n");
   break;
#endif
   // 116 ESTALE           Stale NFS file handle
#ifdef ESTALE
 case ESTALE:
   fprintf(stderr,"116 ESTALE: Stale NFS file handle\n");
   break;
#endif
   // 115 EINPROGRESS      Operation now in progress
#ifdef EINPROGRESS
 case EINPROGRESS:
   fprintf(stderr,"115 EINPROGRESS: Operation now in progress\n");
   break;
#endif
   // 114 EALREADY         Operation already in progress
#ifdef EALREADY
 case EALREADY:
   fprintf(stderr,"114 EALREADY: Operation already in progress\n");
   break;
#endif
   // 113 EHOSTUNREACH     No route to host
#ifdef EHOSTUNREACH
 case EHOSTUNREACH:
   fprintf(stderr,"113 EHOSTUNREACH: No route to host\n");
   break;
#endif
   // 112 EHOSTDOWN        Host is down
#ifdef EHOSTDOWN
 case EHOSTDOWN:
   fprintf(stderr,"112 EHOSTDOWN: Host is down\n");
   break;
#endif
   // 111 ECONNREFUSED     Connection refused
#ifdef ECONNREFUSED
 case ECONNREFUSED:
   fprintf(stderr,"111 ECONNREFUSED: Connection refused\n");
   break;
#endif
   // 110 ETIMEDOUT        Connection timed out
#ifdef ETIMEDOUT
 case ETIMEDOUT:
   fprintf(stderr,"110 ETIMEDOUT: Connection timed out\n");
   break;
#endif
   // 109 ETOOMANYREFS     Too many references: cannot splice
#ifdef ETOOMANYREFS
 case ETOOMANYREFS:
   fprintf(stderr,"109 ETOOMANYREFS: Too many references: cannot splice\n");
   break;
#endif
   // 108 ESHUTDOWN        Cannot send after transport endpoint shutdown
#ifdef ESHUTDOWN
 case ESHUTDOWN:
   fprintf(stderr,"108 ESHUTDOWN: Cannot send after transport endpoint shutdown\n");
   break;
#endif
   // 107 ENOTCONN         Transport endpoint is not connected
#ifdef ENOTCONN
 case ENOTCONN:
   fprintf(stderr,"107 ENOTCONN: Transport endpoint is not connected\n");
   break;
#endif
   // 106 EISCONN          Transport endpoint is already connected
#ifdef EISCONN
 case EISCONN:
   fprintf(stderr,"106 EISCONN: Transport endpoint is already connected\n");
   break;
#endif
   // 105 ENOBUFS          No buffer space available
#ifdef ENOBUFS
 case ENOBUFS:
   fprintf(stderr,"105 ENOBUFS: No buffer space available\n");
   break;
#endif
   // 104 ECONNRESET       Connection reset by peer
#ifdef ECONNRESET
 case ECONNRESET:
   fprintf(stderr,"104 ECONNRESET: Connection reset by peer\n");
   break;
#endif
   // 103 ECONNABORTED     Software caused connection abort
#ifdef ECONNABORTED
 case ECONNABORTED:
   fprintf(stderr,"103 ECONNABORTED: Software caused connection abort\n");
   break;
#endif
   // 102 ENETRESET        Network dropped connection on reset
#ifdef ENETRESET
 case ENETRESET:
   fprintf(stderr,"102 ENETRESET: Network dropped connection on reset\n");
   break;
#endif
   // 101 ENETUNREACH      Network is unreachable
#ifdef ENETUNREACH
 case ENETUNREACH:
   fprintf(stderr,"101 ENETUNREACH: Network is unreachable\n");
   break;
#endif
   // 100 ENETDOWN         Network is down
#ifdef ENETDOWN
 case ENETDOWN:
   fprintf(stderr,"100 ENETDOWN: Network is down\n");
   break;
#endif
   // 99 EADDRNOTAVAIL    Cannot assign requested address
#ifdef ADDRNOTAVAIL
 case EDDRNOTAVAIL:
   fprintf(stderr,"99 EDDRNOTAVAIL: Cannot assign requested address\n");
   break;
#endif
   // 98 EADDRINUSE       Address already in use
#ifdef EADDRINUSE
 case EADDRINUSE:
   fprintf(stderr,"98 EADDRINUSE: Address already in use\n");
   break;
#endif
   // 97 EAFNOSUPPORT     Address family not supported by protocol
#ifdef EAFNOSUPPORT
 case EAFNOSUPPORT:
   fprintf(stderr,"97 EAFNOSUPPORT: Address family not supported by protocol\n");
   break;
#endif
   // 96 EPFNOSUPPORT     Protocol family not supported
#ifdef EPFNOSUPPORT
 case EPFNOSUPPORT:
   fprintf(stderr,"96 EPFNOSUPPORT: Protocol family not supported\n");
   break;
#endif
   // 95 EOPNOTSUPP       Operation not supported
#ifdef EOPNOTSUPP
 case EOPNOTSUPP:
   fprintf(stderr,"95 EOPNOTSUPP: Operation not supported\n");
   break;
#endif
   // 94 ESOCKTNOSUPPORT  Socket type not supported
#ifdef ESOCKTNOSUPPORT
 case ESOCKTNOSUPPORT:
   fprintf(stderr,"94 ESOCKTNOSUPPORT: Socket type not supported\n");
   break;
#endif
   // 93 EPROTONOSUPPORT  Protocol not supported
#ifdef EPROTONOSUPPORT
 case EPROTONOSUPPORT:
   fprintf(stderr,"93 EPROTONOSUPPORT: Protocol not supported\n");
   break;
#endif
   // 92 ENOPROTOOPT      Protocol not available
#ifdef ENOPROTOOPT
 case ENOPROTOOPT:
   fprintf(stderr,"92 ENOPROTOOPT: Protocol not available\n");
   break;
#endif
   // 91 EPROTOTYPE       Protocol wrong type for socket
#ifdef EPROTOTYPE
 case EPROTOTYPE:
   fprintf(stderr,"91 EPROTOTYPE: Protocol wrong type for socket\n");
   break;
#endif
   // 90 EMSGSIZE         Message too long
#ifdef EMSGSIZE
 case EMSGSIZE:
   fprintf(stderr,"90 EMSGSIZE: Message too long\n");
   break;
#endif
   // 89 EDESTADDRREQ     Destination address required
#ifdef EDESTADDRREQ
 case EDESTADDRREQ:
   fprintf(stderr,"89 EDESTADDRREQ: Destination address required\n");
   break;
#endif
   // 88 ENOTSOCK         Socket operation on non-socket
#ifdef ENOTSOCK
 case ENOTSOCK:
   fprintf(stderr,"88 ENOTSOCK: Socket operation on non-socket\n");
   break;
#endif
   // 87 EUSERS           Too many users
#ifdef EUSERS
 case EUSERS:
   fprintf(stderr,"87 EUSERS: Too many users\n");
   break;
#endif
   // 86 ESTRPIPE         Streams pipe error
#ifdef ESTRPIPE
 case ESTRPIPE:
   fprintf(stderr,"86 ESTRPIPE: Streams pipe error\n");
   break;
#endif
   // 85 ERESTART         Interrupted system call should be restarted
#ifdef ERESTART
 case ERESTART:
   fprintf(stderr,"85 ERESTART: Interrupted system call should be restarted\n");
   break;
#endif
   // 84 EILSEQ           Invalid or incomplete multibyte or wide character
#ifdef EILSEQ
 case EILSEQ:
   fprintf(stderr,"84 EILSEQ: Invalid or incomplete multibyte or wide character\n");
   break;
#endif
   // 83 ELIBEXEC         Cannot exec a shared library directly
#ifdef ELIBEXEC
 case ELIBEXEC:
   fprintf(stderr,"83 ELIBEXEC: Cannot exec a shared library directly\n");
   break;
#endif
   // 82 ELIBMAX          Attempting to link in too many shared libraries
#ifdef ELIBMAX
 case ELIBMAX:
   fprintf(stderr,"82 ELIBMAX: Attempting to link in too many shared libraries\n");
   break;
#endif
   // 81 ELIBSCN          .lib section in a.out corrupted
#ifdef ELIBSCN
 case ELIBSCN:
   fprintf(stderr,"81 ELIBSCN: .lib section in a.out corrupted\n");
   break;
#endif
   // 80 ELIBBAD          Accessing a corrupted shared library
#ifdef ELIBBAD
 case ELIBBAD:
   fprintf(stderr,"80 ELIBBAD: Accessing a corrupted shared library\n");
   break;
#endif
   // 79 ELIBACC          Can not access a needed shared library
#ifdef ELIBACC
 case ELIBACC:
   fprintf(stderr,"79 ELIBACC: Can not access a needed shared library\n");
   break;
#endif
   // 78 EREMCHG          Remote address changed
#ifdef EREMCHG
 case EREMCHG:
   fprintf(stderr,"78 EREMCHG: Remote address changed\n");
   break;
#endif
   // 77 EBADFD           File descriptor in bad state
#ifdef EBADFD
 case EBADFD:
   fprintf(stderr,"77 EBADFD: File descriptor in bad state\n");
   break;
#endif
   // 76 ENOTUNIQ         Name not unique on network
#ifdef ENOTUNIQ
 case ENOTUNIQ:
   fprintf(stderr,"76 ENOTUNIQ: Name not unique on network\n");
   break;
#endif
   // 75 EOVERFLOW        Value too large for defined data type
#ifdef EOVERFLOW
 case EOVERFLOW:
   fprintf(stderr,"75 EOVERFLOW: Value too large for defined data type\n");
   break;
#endif
   // 74 EBADMSG          Bad message
#ifdef EBADMSG
 case EBADMSG:
   fprintf(stderr,"74 EBADMSG: Bad message\n");
   break;
#endif
   // 73 EDOTDOT          RFS specific error
#ifdef EDOTDOT
 case EDOTDOT:
   fprintf(stderr,"73 EDOTDOT: RFS specific error\n");
   break;
#endif
   // 72 EMULTIHOP        Multihop attempted
#ifdef EMULTIHOP
 case EMULTIHOP:
   fprintf(stderr,"72 EMULTIHOP: Multihop attempted\n");
   break;
#endif
   // 71 EPROTO           Protocol error
#ifdef EPROTO
 case EPROTO:
   fprintf(stderr,"71 EPROTO: Protocol error\n");
   break;
#endif
   // 70 ECOMM            Communication error on send
#ifdef ECOMM
 case ECOMM:
   fprintf(stderr,"70 ECOMM: Communication error on send\n");
   break;
#endif
   // 69 ESRMNT           Srmount error
#ifdef ESRMNT
 case ESRMNT:
   fprintf(stderr,"69 ESRMNT: Srmount error\n");
   break;
#endif
   // 68 EADV             Advertise error
#ifdef EADV
 case EADV:
   fprintf(stderr,"68 EADV: Advertise error\n");
   break;
#endif
   // 67 ENOLINK          Link has been severed
#ifdef ENOLINK
 case ENOLINK:
   fprintf(stderr,"67 ENOLINK: Link has been severed\n");
   break;
#endif
   // 66 EREMOTE          Object is remote
#ifdef EREMOTE
 case EREMOTE:
   fprintf(stderr,"66 EREMOTE: Object is remote\n");
   break;
#endif
   // 65 ENOPKG           Package not installed
#ifdef ENOPKG
 case ENOPKG:
   fprintf(stderr,"65 ENOPKG: Package not installed\n");
   break;
#endif
   // 64 ENONET           Machine is not on the network
#ifdef ENONET
 case ENONET:
   fprintf(stderr,"64 ENONET: Machine is not on the network\n");
   break;
#endif
   // 63 ENOSR            Out of streams resources
#ifdef ENOSR
 case ENOSR:
   fprintf(stderr,"63 ENOSR: Out of streams resources\n");
   break;
#endif
   // 62 ETIME            Timer expired
#ifdef ETIME
 case ETIME:
   fprintf(stderr,"62 ETIME: Timer expired\n");
   break;
#endif
   // 61 ENODATA          No data available
#ifdef ENODATA
 case ENODATA:
   fprintf(stderr,"61 ENODATA: No data available\n");
   break;
#endif
   // 60 ENOSTR           Device not a stream
#ifdef ENOSTR
 case ENOSTR:
   fprintf(stderr,"60 ENOSTR: Device not a stream\n");
   break;
#endif
   // 59 EBFONT           Bad font file format
#ifdef EBFONT
 case EBFONT:
   fprintf(stderr,"59 EBFONT: Bad font file format\n");
   break;
#endif
   // 57 EBADSLT          Invalid slot
#ifdef EBADSLT
 case EBADSLT:
   fprintf(stderr,"57 EBADSLT: Invalid slot\n");
   break;
#endif
   // 56 EBADRQC          Invalid request code
#ifdef EBADRQC
 case EBADRQC:
   fprintf(stderr,"56 EBADRQC: Invalid request code\n");
   break;
#endif
   // 55 ENOANO           No anode
#ifdef ENOANO
 case ENOANO:
   fprintf(stderr,"55 ENOANO: No anode\n");
   break;
#endif
   // 54 EXFULL           Exchange full
#ifdef EXFULL
 case EXFULL:
   fprintf(stderr,"54 EXFULL: Exchange full\n");
   break;
#endif
   // 53 EBADR            Invalid request descriptor
#ifdef EBADR
 case EBADR:
   fprintf(stderr,"53 EBADR: Invalid request descriptor\n");
   break;
#endif
   // 52 EBADE            Invalid exchange
#ifdef EBADE
 case EBADE:
   fprintf(stderr,"52 EBADE: Invalid exchange\n");
   break;
#endif
   // 51 EL2HLT           Level 2 halted
#ifdef EL2HLT
 case EL2HLT:
   fprintf(stderr,"51 EL2HLT: Level 2 halted\n");
   break;
#endif
   // 50 ENOCSI           No CSI structure available
#ifdef ENOCSI
 case ENOCSI:
   fprintf(stderr,"50 ENOCSI: No CSI structure available\n");
   break;
#endif
   // 49 EUNATCH          Protocol driver not attached
#ifdef EUNATCH
 case EUNATCH:
   fprintf(stderr,"49 EUNATCH: Protocol driver not attached\n");
   break;
#endif
   // 48 ELNRNG           Link number out of range
#ifdef ELNRNG
 case ELNRNG:
   fprintf(stderr,"48 ELNRNG: Link number out of range\n");
   break;
#endif
   // 47 EL3RST           Level 3 reset
#ifdef EL3RST
 case EL3RST:
   fprintf(stderr,"47 EL3RST: Level 3 reset\n");
   break;
#endif
   // 46 EL3HLT           Level 3 halted
#ifdef EL3HLT
 case EL3HLT:
   fprintf(stderr,"46 EL3HLT: Level 3 halted\n");
   break;
#endif
   // 45 EL2NSYNC         Level 2 not synchronized
#ifdef EL2NSYNC
 case EL2NSYNC:
   fprintf(stderr,"45 EL2NSYNC: Level 2 not synchronized\n");
   break;
#endif
   // 44 ECHRNG           Channel number out of range
#ifdef ECHRNG
 case ECHRNG:
   fprintf(stderr,"44 ECHRNG: Channel number out of range\n");
   break;
#endif
   // 43 EIDRM            Identifier removed
#ifdef EIDRM
 case EIDRM:
   fprintf(stderr,"43 EIDRM: Identifier removed\n");
   break;
#endif
   // 42 ENOMSG           No message of desired type
#ifdef ENOMSG
 case ENOMSG:
   fprintf(stderr,"42 ENOMSG: No message of desired type\n");
   break;
#endif
   // 40 ELOOP            Too many levels of symbolic links
#ifdef ELOOP
 case ELOOP:
   fprintf(stderr,"40 ELOOP: Too many levels of symbolic links\n");
   break;
#endif
   // 39 ENOTEMPTY        Directory not empty
#ifdef ENOTEMPTY
 case ENOTEMPTY:
   fprintf(stderr,"39 ENOTEMPTY: Directory not empty\n");
   break;
#endif
   // 38 ENOSYS           Function not implemented
#ifdef ENOSYS
 case ENOSYS:
   fprintf(stderr,"38 ENOSYS: Function not implemented\n");
   break;
#endif
   // 37 ENOLCK           No locks available
#ifdef ENOLCK
 case ENOLCK:
   fprintf(stderr,"37 ENOLCK: No locks available\n");
   break;
#endif
   // 36 ENAMETOOLONG     File name too long
#ifdef ENAMETOOLONG
 case ENAMETOOLONG:
   fprintf(stderr,"36 ENAMETOOLONG: File name too long\n");
   break;
#endif
   // 35 EDEADLK          Resource deadlock avoided
#ifdef EDEADLK
 case EDEADLK:
   fprintf(stderr,"35 EDEADLK: Resource deadlock avoided\n");
   break;
#endif
   // 34 ERANGE           Numerical result out of range
#ifdef ERANGE
 case ERANGE:
   fprintf(stderr,"34 ERANGE: Numerical result out of range\n");
   break;
#endif
   // 33 EDOM             Numerical argument out of domain
#ifdef EDOM
 case EDOM:
   fprintf(stderr,"33 EDOM: Numerical argument out of domain\n");
   break;
#endif
   // 32 EPIPE            Broken pipe
#ifdef EPIPE
 case EPIPE:
   fprintf(stderr,"32 EPIPE: Broken pipe\n");
   break;
#endif
   // 31 EMLINK           Too many links
#ifdef EMLINK
 case EMLINK:
   fprintf(stderr,"31 EMLINK: Too many links\n");
   break;
#endif
   // 30 EROFS            Read-only file system
#ifdef EROFS
 case EROFS:
   fprintf(stderr,"30 EROFS: Read-only file system\n");
   break;
#endif
   // 29 ESPIPE           Illegal seek
#ifdef ESPIPE
 case ESPIPE:
   fprintf(stderr,"29 ESPIPE: Illegal seek\n");
   break;
#endif
   // 28 ENOSPC           No space left on device
#ifdef ENOSPC
 case ENOSPC:
   fprintf(stderr,"28 ENOSPC: No space left on device\n");
   break;
#endif
   // 27 EFBIG            File too large
#ifdef EFBIG
 case EFBIG:
   fprintf(stderr,"27 EFBIG: File too large\n");
   break;
#endif
   // 26 ETXTBSY          Text file busy
#ifdef ETXTBSY
 case ETXTBSY:
   fprintf(stderr,"26 ETXTBSY: Text file busy\n");
   break;
#endif
   // 25 ENOTTY           Inappropriate ioctl for device
#ifdef ENOTTY
 case ENOTTY:
   fprintf(stderr,"25 ENOTTY: Inappropriate ioctl for device\n");
   break;
#endif
   // 24 EMFILE           Too many open files
#ifdef EMFILE
 case EMFILE:
   fprintf(stderr,"24 EMFILE: Too many open files\n");
   break;
#endif
   // 23 ENFILE           Too many open files in system
#ifdef ENFILE
 case ENFILE:
   fprintf(stderr,"23 ENFILE: Too many open files in system\n");
   break;
#endif
   // 22 EINVAL           Invalid argument
#ifdef EINVAL
 case EINVAL:
   fprintf(stderr,"22 EINVAL: Invalid argument\n");
   break;
#endif
   // 21 EISDIR           Is a directory
#ifdef EISDIR
 case EISDIR:
   fprintf(stderr,"21 EISDIR: Is a directory\n");
   break;
#endif
   // 20 ENOTDIR          Not a directory
#ifdef ENOTDIR
 case ENOTDIR:
   fprintf(stderr,"20 ENOTDIR: Not a directory\n");
   break;
#endif
   // 19 ENODEV           No such device
#ifdef ENODEV
 case ENODEV:
   fprintf(stderr,"19 ENODEV: No such device\n");
   break;
#endif
   // 18 EXDEV            Invalid cross-device link
#ifdef EXDEV
 case EXDEV:
   fprintf(stderr,"18 EXDEV: Invalid cross-device link\n");
   break;
#endif
   // 17 EEXIST           File exists
#ifdef EEXIST
 case EEXIST:
   fprintf(stderr,"17 EEXIST: File exists\n");
   break;
#endif
   // 16 EBUSY            Device or resource busy
#ifdef EBUSY
 case EBUSY:
   fprintf(stderr,"16 EBUSY: Device or resource busy\n");
   break;
#endif
   // 15 ENOTBLK          Block device required
#ifdef ENOTBLK
 case ENOTBLK:
   fprintf(stderr,"15 ENOTBLK: Block device required\n");
   break;
#endif
   // 14 EFAULT           Bad address
#ifdef EFAULT
 case EFAULT:
   fprintf(stderr,"14 EFAULT: Bad address\n");
   break;
#endif
   // 13 EACCES           Permission denied
#ifdef EACCES
 case EACCES:
   fprintf(stderr,"13 EACCES: Permission denied\n");
   break;
#endif
   // 12 ENOMEM           Cannot allocate memory
#ifdef ENOMEM
 case ENOMEM:
   fprintf(stderr,"12 ENOMEM: Cannot allocate memory\n");
   break;
#endif
   // 11 EAGAIN           Resource temporarily unavailable
#ifdef EAGAIN
 case EAGAIN:
   fprintf(stderr,"11 EAGAIN: Resource temporarily unavailable\n");
   break;
#endif
   // 10 ECHILD           No child processes
#ifdef ECHILD
 case ECHILD:
   fprintf(stderr,"10 ECHILD: No child processes\n");
   break;
#endif
   // 9 EBADF            Bad file descriptor
#ifdef EBADF
 case EBADF:
   fprintf(stderr,"9 EBADF: Bad file descriptor\n");
   break;
#endif
   // 8 ENOEXEC          Exec format error
#ifdef ENOEXEC
 case ENOEXEC:
   fprintf(stderr,"8 ENOEXEC: Exec format error\n");
   break;
#endif
   // 7 E2BIG            Argument list too long
#ifdef E2BIG
 case E2BIG:
   fprintf(stderr,"7 E2BIG: Argument list too long\n");
   break;
#endif
   // 6 ENXIO            No such device or address
#ifdef ENXIO
 case ENXIO:
   fprintf(stderr,"6 ENXIO: No such device or address\n");
   break;
#endif
   // 5 EIO              Input/output error
#ifdef EIO
 case EIO:
   fprintf(stderr,"5 EIO: Input/output error\n");
   break;
#endif
   // 4 EINTR            Interrupted system call
#ifdef EINTR
 case EINTR:
   fprintf(stderr,"4 EINTR: Interrupted system call\n");
   break;
#endif
   // 3 ESRCH            No such process
#ifdef ESRCH
 case ESRCH:
   fprintf(stderr,"3 ESRCH: No such process\n");
   break;
#endif
   // 2 ENOENT           No such file or directory
#ifdef ENOENT
 case ENOENT:
   fprintf(stderr,"2 ENOENT: No such file or directory\n");
   break;
#endif
   // 1 EPERM            Operation not permitted
#ifdef EPERM
 case EPERM:
   fprintf(stderr,"1 EPERM: Operation not permitted\n");
   break;
#endif

  default:
    fprintf(stderr,"Unknown errno\n");
    exit(EX_USAGE);
  }
}
