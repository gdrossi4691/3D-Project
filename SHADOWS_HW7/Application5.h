// Application4.h: interface for the Application4 class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_APPLICATION5_H__43A7FA9C_6CD6_4A79_9567_2354BFEFAFFB__INCLUDED_)
#define AFX_APPLICATION5_H__43A7FA9C_6CD6_4A79_9567_2354BFEFAFFB__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

// for HW6
#define	AAKERNEL_SIZE	6

#include "Application.h"


class Application5 : public Application  
{
private:
	GzDisplay*  m_pDisplays[AAKERNEL_SIZE];
public:
	Application5();
	virtual ~Application5();
	
	int	Initialize();
	int LoadModel(GzCoord* vertexLists, GzCoord* normalLists, GzTextureIndex* uvLists);
	int LoadObjModel(Model** triangles);
	virtual int Render(); 
	int Clean();
};

#endif // !defined(AFX_APPLICATION5_H__43A7FA9C_6CD6_4A79_9567_2354BFEFAFFB__INCLUDED_)
