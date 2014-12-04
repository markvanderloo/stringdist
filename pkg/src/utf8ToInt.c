/*  stringdist - a C library of string distance algorithms with an interface to R.
 *  Copyright (C) 2013  Mark van der Loo
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>. 
 *
 *  You can contact the author at: mark _dot_ vanderloo _at_ gmail _dot_ com
 */



/* This function is gratefully copied from the R core distribution.
 * It is replicated here to facilitate multicore processing through
 * openmp (as it is not part of R's API).
 *
 */
static int mbrtoint(int *w, const char *s)
{
    unsigned int byte;
    byte = *((unsigned char *)s);

    if (byte == 0) {
  *w = 0;
  return 0;
    } else if (byte < 0xC0) {
  *w = (int) byte;
  return 1;
    } else if (byte < 0xE0) {
  if (!s[1]) return -2;
  if ((s[1] & 0xC0) == 0x80) {
      *w = (int) (((byte & 0x1F) << 6) | (s[1] & 0x3F));
      return 2;
  } else return -1;
    } else if (byte < 0xF0) {
  if (!s[1] || !s[2]) return -2;
  if (((s[1] & 0xC0) == 0x80) && ((s[2] & 0xC0) == 0x80)) {
      *w = (int) (((byte & 0x0F) << 12)
      | ((s[1] & 0x3F) << 6) | (s[2] & 0x3F));
      byte = *w;
      if (byte >= 0xD800 && byte <= 0xDFFF) return -1; /* surrogate */
      if (byte == 0xFFFE || byte == 0xFFFF) return -1;
      return 3;
  } else return -1;
    } else if (byte < 0xF8) {
  if (!s[1] || !s[2] || !s[3]) return -2;
  if (((s[1] & 0xC0) == 0x80)
      && ((s[2] & 0xC0) == 0x80)
      && ((s[3] & 0xC0) == 0x80)) {
      *w = (int) (((byte & 0x07) << 18)
      | ((s[1] & 0x3F) << 12)
      | ((s[2] & 0x3F) << 6)
      | (s[3] & 0x3F));
      byte = *w;
      return 4;
  } else return -1;
    } else if (byte < 0xFC) {
  if (!s[1] || !s[2] || !s[3] || !s[4]) return -2;
  if (((s[1] & 0xC0) == 0x80)
      && ((s[2] & 0xC0) == 0x80)
      && ((s[3] & 0xC0) == 0x80)
      && ((s[4] & 0xC0) == 0x80)) {
      *w = (int) (((byte & 0x03) << 24)
      | ((s[1] & 0x3F) << 18)
      | ((s[2] & 0x3F) << 12)
      | ((s[3] & 0x3F) << 6)
      | (s[4] & 0x3F));
      byte = *w;
      return 5;
  } else return -1;
    } else {
  if (!s[1] || !s[2] || !s[3] || !s[4] || !s[5]) return -2;
  if (((s[1] & 0xC0) == 0x80)
      && ((s[2] & 0xC0) == 0x80)
      && ((s[3] & 0xC0) == 0x80)
      && ((s[4] & 0xC0) == 0x80)
      && ((s[5] & 0xC0) == 0x80)) {
      *w = (int) (((byte & 0x01) << 30)
      | ((s[1] & 0x3F) << 24)
      | ((s[2] & 0x3F) << 18)
      | ((s[3] & 0x3F) << 12)
      | ((s[5] & 0x3F) << 6)
      | (s[5] & 0x3F));
      byte = *w;
      return 6;
  } else return -1;
    }
    /* return -2; not reached */
}

