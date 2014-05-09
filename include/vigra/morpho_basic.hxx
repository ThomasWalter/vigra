
// Author(s): Thomas Walter
// $Date$
// $Rev$
// $URL$

#ifndef MORPHO_BASIC_HXX_
#define MORPHO_BASIC_HXX_


#include "morpho_utilities.hxx"


namespace vigra{
namespace morpho{

template<class Iterator1, class Accessor1, class Iterator2, class Accessor2, class SElement, class Functor>
void morphoBasicSEOperation( Iterator1 srcUpperLeft, Iterator1 srcLowerRight, Accessor1 srca,
                    Iterator2 destUpperLeft, Accessor2 desta,
                    SElement & se,
                    Functor f) {
  typename Accessor1::value_type neutralElement = f.neutralValue;
  typename Accessor1::value_type localMax;

  // border treatment
  // offsets correspond to the maximal extension of the SE.
  Diff2D minOffset = se.minOffset();
  Diff2D maxOffset = se.maxOffset();

  const Iterator1 upperLeftCorner = srcUpperLeft;
  const Iterator1 lowerRightCorner = srcLowerRight;
  const Iterator1 upperLeftCritical = srcUpperLeft - minOffset;
  const Iterator1 lowerRightCritical = srcLowerRight - maxOffset;

  for(; srcUpperLeft.y < upperLeftCritical.y; ++srcUpperLeft.y, ++destUpperLeft.y)
  {
      Iterator1 scurrent = srcUpperLeft;
      Iterator2 dcurrent = destUpperLeft;

      for(; scurrent.x < srcLowerRight.x; ++scurrent.x, ++dcurrent.x)
      {
          localMax = neutralElement;
          for(typename SElement::ITERATORTYPE iter = se.begin();
              iter != se.end();
              ++iter)
          {
              if(    ( (scurrent + *iter).y >= upperLeftCorner.y)
                  && ( (scurrent + *iter).x >= upperLeftCorner.x)
                  && ( (scurrent + *iter).x < lowerRightCorner.x))
                  localMax = f(localMax, srca(scurrent + *iter));
          }
          desta.set(localMax, dcurrent);
      } // end of x loop
  } // end for the first y-loop.

  for(; srcUpperLeft.y < lowerRightCritical.y; ++srcUpperLeft.y, ++destUpperLeft.y)
  {
      Iterator1 scurrent = srcUpperLeft;
      Iterator2 dcurrent = destUpperLeft;

      // x-loop: the left side
      for(; scurrent.x < upperLeftCritical.x; ++scurrent.x, ++dcurrent.x)
      {
          localMax = neutralElement;
          for(typename SElement::ITERATORTYPE iter = se.begin();
              iter != se.end();
              ++iter)
          {
              if( (scurrent + *iter).x >= upperLeftCorner.x )
                  localMax = f(localMax, srca(scurrent + *iter));
          }
          desta.set(localMax, dcurrent);
      } // end of x loop (left)

      for(; scurrent.x < lowerRightCritical.x; ++scurrent.x, ++dcurrent.x)
      {
          localMax = neutralElement;
          for(typename SElement::ITERATORTYPE iter = se.begin();
              iter != se.end();
              ++iter)
          {
              localMax = f(localMax, srca(scurrent + *iter));
          }
          desta.set(localMax, dcurrent);
      } // end of the middle x loop

      // the right side
      for(; scurrent.x < srcLowerRight.x; ++scurrent.x, ++dcurrent.x)
      {
          localMax = neutralElement;
          for(typename SElement::ITERATORTYPE iter = se.begin();
              iter != se.end();
              ++iter)
          {
              if( (scurrent + *iter).x < lowerRightCorner.x)
                  localMax = f(localMax, srca(scurrent + *iter));
          }
          desta.set(localMax, dcurrent);
      } // end of the right x loop
  } // end of y loop (middle)

  // y-loop: lower
  for(; srcUpperLeft.y < srcLowerRight.y; ++srcUpperLeft.y, ++destUpperLeft.y)
  {
      Iterator1 scurrent = srcUpperLeft;
      Iterator2 dcurrent = destUpperLeft;

      for(; scurrent.x < srcLowerRight.x; ++scurrent.x, ++dcurrent.x)
      {
          localMax = neutralElement;
          for(typename SElement::ITERATORTYPE iter = se.begin();
              iter != se.end();
              ++iter)
          {
              if(    ( (scurrent + *iter).y < lowerRightCorner.y)
                  && ( (scurrent + *iter).x < lowerRightCorner.x)
                  && ( (scurrent + *iter).x >= upperLeftCorner.x) )
                  localMax = f(localMax, srca(scurrent + *iter));
          }
          desta.set(localMax, dcurrent);
      } // end of x loop
  } // end for the lower y-loop.
} // end of morphoBasicSEOperation

/////////////////////////////////////////////////////////////////////////
// EROSION AND DILATION
/////////////////////////////////////////////////////////////////////////

// Morphological dilation
template<class Iterator1, class Accessor1, class Iterator2, class Accessor2, class SElement>
void morphoDilation(vigra::triple<Iterator1, Iterator1, Accessor1> src,
              vigra::triple<Iterator2, Iterator2, Accessor2> dest,
              SElement se)
{
    vigra::BasicImage<typename Accessor2::value_type> temp(dest.second - dest.first);

	if(se.size > 0)
    morphoBasicSEOperation(src.first, src.second, src.third,
                     dest.first, dest.third,
                      se,
                     MaxFunctor<typename Accessor1::value_type>());

  // a morphological dilation with se of size n
  // corresponds to n morphological dilations with size 1.
  for(int i = 1; i < se.size; i++)
  {
	  if(i%2 == 0)
		  morphoBasicSEOperation(temp.upperLeft(), temp.lowerRight(), temp.accessor(),
		                dest.first, dest.third,
		                se,
		                MaxFunctor<typename Accessor1::value_type>());
	  else
		  morphoBasicSEOperation(dest.first, dest.second, dest.third,
		                temp.upperLeft(), temp.accessor(),
		                se,
		                MaxFunctor<typename Accessor1::value_type>());
  }

       if(se.size%2 == 0)
    	   vigra::copyImage(temp.upperLeft(), temp.lowerRight(), temp.accessor(),
    	                    dest.first, dest.third);
} // end of dilation

// Morphological erosion
template<class Iterator1, class Accessor1, class Iterator2, class Accessor2,class SElement>
void morphoErosion(vigra::triple<Iterator1, Iterator1, Accessor1> src,
         vigra::triple<Iterator2, Iterator2, Accessor2> dest,
         SElement se)
{

	vigra::BasicImage<typename Accessor2::value_type> temp(dest.second - dest.first);

	if(se.size > 0)
		morphoBasicSEOperation(src.first, src.second, src.third,
		              dest.first, dest.third,
		              se,
		              MinFunctor<typename Accessor1::value_type>());

  // a morphological erosion with se of size n
  // corresponds to n morphological erosions with size 1.
  for(int i = 1; i < se.size; i++)
  {
    if(i%2 == 0)
    	morphoBasicSEOperation(temp.upperLeft(), temp.lowerRight(), temp.accessor(),
    	              dest.first, dest.third,
    	              se,
    	              MinFunctor<typename Accessor1::value_type>());
    else
    	morphoBasicSEOperation(dest.first, dest.second, dest.third,
    	              temp.upperLeft(), temp.accessor(),
    	              se,
    	              MinFunctor<typename Accessor1::value_type>());
  }

  if(se.size%2 == 0)
	  vigra::copyImage(temp.upperLeft(), temp.lowerRight(), temp.accessor(),
			  dest.first, dest.third);

} // end of erosion


/////////////////////////////////////////////////////////////////////////
// OPENING AND CLOSING
/////////////////////////////////////////////////////////////////////////

// Morphological opening
template<class Iterator1, class Accessor1,
     class Iterator2, class Accessor2,
     class SElement>
void morphoOpening(vigra::triple<Iterator1, Iterator1, Accessor1> src,
            vigra::triple<Iterator2, Iterator2, Accessor2> dest, SElement se)
{
    vigra::BasicImage<typename Accessor1::value_type> temp(src.second - src.first);
    morphoErosion(src, vigra::destImageRange(temp), se);
    se.transpose();
    morphoDilation(vigra::srcImageRange(temp), dest, se);
    se.transpose();
} // end of opening

// Morphological closing
template<class Iterator1, class Accessor1,
         class Iterator2, class Accessor2,
         class SElement>
void morphoClosing(vigra::triple<Iterator1, Iterator1, Accessor1> src,
             vigra::triple<Iterator2, Iterator2, Accessor2> dest, SElement se)
{
    vigra::BasicImage<typename Accessor1::value_type> temp(src.second - src.first);
    morphoDilation(src, vigra::destImageRange(temp), se);
    se.transpose();
    morphoErosion(vigra::srcImageRange(temp), dest, se);
    se.transpose();
} // end of closing

// Morphological gradient
template<class Iterator1, class Accessor1,
         class Iterator2, class Accessor2,
         class SElement>
void morphoGradient(vigra::triple<Iterator1, Iterator1, Accessor1> src,
                      vigra::triple<Iterator2, Iterator2, Accessor2> dest,
                      SElement se,
                      typename Accessor2::value_type markerVal = 255)
{
    typedef typename Accessor1::value_type INTYPE;
    typedef typename Accessor2::value_type OUTTYPE;

    vigra::BasicImage<INTYPE> dil(src.second - src.first);
    vigra::BasicImage<INTYPE> ero(src.second - src.first);
    morphoDilation(src, vigra::destImageRange(dil), se);
    morphoErosion(src, vigra::destImageRange(ero), se);
    vigra::combineTwoImages(srcImageRange(dil), srcImage(ero),
                            destIter(dest.first, dest.third),
                            std::minus<INTYPE>() );

} // end of morphogradient

// External gradient
template<class Iterator1, class Accessor1,
     class Iterator2, class Accessor2,
     class SElement>
void morphoExternalGradient(vigra::triple<Iterator1, Iterator1, Accessor1> src,
                        vigra::triple<Iterator2, Iterator2, Accessor2> dest,
                        SElement se,
                        typename Accessor2::value_type markerVal = 255)
{
    typedef typename Accessor1::value_type INTYPE;
    typedef typename Accessor2::value_type OUTTYPE;

    vigra::BasicImage<INTYPE> dil(src.second - src.first);
    morphoDilation(src, vigra::destImageRange(dil), se);
    vigra::combineTwoImages(srcImageRange(dil),
                            srcIter(src.first, src.third),
                            destIter(dest.first, dest.third),
                            std::minus<INTYPE>() );

}

// Internal gradient
template<class Iterator1, class Accessor1,
         class Iterator2, class Accessor2,
         class SElement>
void morphoInternalGradient(vigra::triple<Iterator1, Iterator1, Accessor1> src,
                        vigra::triple<Iterator2, Iterator2, Accessor2> dest,
                        SElement se,
                        typename Accessor2::value_type markerVal = 255)
{
    typedef typename Accessor1::value_type INTYPE;
    typedef typename Accessor2::value_type OUTTYPE;

    vigra::BasicImage<INTYPE> ero(src.second - src.first);
    morphoErosion(src, vigra::destImageRange(ero), se);
    vigra::combineTwoImages(srcIterRange(src.first, src.second, src.third),
                            srcImage(ero),
                            destIter(dest.first, dest.third),
                            std::minus<INTYPE>() );

}


template<class Image1, class Image2, class SElement>
void morphoInternalGradient(const Image1 & imin, Image2 & imout, SElement se)
{
  morphoInternalGradient(srcImageRange(imin), destImageRange(imout), se);
}

template<class Image1, class Image2, class SElement>
void morphoExternalGradient(const Image1 & imin, Image2 & imout, SElement se)
{
  morphoExternalGradient(srcImageRange(imin), destImageRange(imout), se);
}

template<class Image1, class Image2, class SElement>
void morphoGradient(const Image1 & imin, Image2 & imout, SElement se)
{
  morphoGradient(srcImageRange(imin), destImageRange(imout), se);
}

template<class Image1, class Image2, class SElement>
void morphoErosion(const Image1 & imin, Image2 & imout, SElement se)
{
  morphoErosion(srcImageRange(imin), destImageRange(imout), se);
}

template<class Image1, class Image2, class SElement>
void morphoDilation(const Image1 & imin, Image2 & imout, SElement se)
{
  morphoDilation(srcImageRange(imin), destImageRange(imout), se);
}

// Open and close
template<class Image1, class Image2, class SElement>
void morphoOpening(const Image1 & imin, Image2 & imout, SElement se)
{
  morphoOpening(srcImageRange(imin), destImageRange(imout), se);
}

template<class Image1, class Image2, class SElement>
void morphoClosing(const Image1 & imin, Image2 & imout, SElement se)
{
  morphoClosing(srcImageRange(imin), destImageRange(imout), se);
}


}; /* end of namespace morpho */
}; /* end of namespace vigra */

#endif /*BASIC_MORPHO_HXX_*/
