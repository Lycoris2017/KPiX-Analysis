//-----------------------------------------------------------------------------
// File          : analysis.h
// Author        : Mengqing Wu (orig. Uwe, Kraemer & Ryan,Herbst) <mengqing.wu@desy.de>
// Created       : 08/16/2018
// Project       : KPiX Analysis
//-----------------------------------------------------------------------------
// Description :
// General analysis of KPiX data.
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
// Modification history :
// 16/08/2018: refactorize kpix analysis code (orig. Ryan Herbst, then modified by Uwe Kraemer)
//-----------------------------------------------------------------------------

#include <DataRead.h>


bool gtx_ltz ( int x, int y, int z) // simple function asking if x<y<z
{
	return ((x <= y) && (y <= z));
}

void addDoubleToXml ( ofstream *xml, uint indent, string variable, Double_t value ) { //xml output function
	uint x;
	
	if ( !(value != value) ) {
		for (x=0; x < indent; x++) *xml << " ";
		*xml << "<" << variable << ">";
		*xml << value;
		*xml << "</" << variable << ">";
		*xml << endl;
	}
}

void addStringToXml ( ofstream *xml, uint indent, string variable, string value ) { //xml output function
	uint x;
	
	for (x=0; x < indent; x++) *xml << " ";
	*xml << "<" << variable << ">";
	*xml << value;
	*xml << "</" << variable << ">";
	*xml << endl;
}

// Class definition
namespace lycoris
{
  class kpixAna {
    
  protected:
    DataRead               dataRead;  //kpix event classes used for analysis of binary date
    
  public:
    
    
  };
  

}// namespace lycoris end
