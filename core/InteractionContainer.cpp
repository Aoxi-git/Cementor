// 2008 © Sergei Dorofeenko <sega@users.berlios.de>
// 2009,2010 © Václav Šmilauer <eudoxos@arcig.cz>

#include "InteractionContainer.hpp"
#include "Scene.hpp"

#ifdef YADE_OPENMP
	#include<omp.h>
#endif
CREATE_LOGGER(InteractionContainer);
// begin internal functions

bool InteractionContainer::insert(const shared_ptr<Interaction>& i){
	assert(bodies);
	boost::mutex::scoped_lock lock(drawloopmutex);

	Body::id_t id1=i->getId1();
	Body::id_t id2=i->getId2();

	if (id1>id2) swap(id1,id2);

	assert((Body::id_t)bodies->size()>id1); // the bodies must exist already
	assert((Body::id_t)bodies->size()>id2);

	const shared_ptr<Body>& b1=(*bodies)[id1];
	const shared_ptr<Body>& b2=(*bodies)[id2];

	if(!b1->intrs.insert(Body::MapId2IntrT::value_type(id2,i)).second) return false; // already exists
	if(!b2->intrs.insert(Body::MapId2IntrT::value_type(id1,i)).second) return false;

	linIntrs.resize(++currSize); // currSize updated
	linIntrs[currSize-1]=i; // assign last element
	i->linIx=currSize-1; // store the index back-reference in the interaction (so that it knows how to erase/move itself)

	const shared_ptr<Scene>& scene=Omega::instance().getScene();
	i->iterBorn=scene->iter;

	return true;
}




bool InteractionContainer::insertInteractionMPI(shared_ptr<Interaction>& i){
// 	if (! i->isReal()){return false; }
	bool res=true;
	const shared_ptr<Interaction>& iOld =  find(Body::id_t(i->id1),Body::id_t(i->id2));
	if (not iOld) res = insert(i);// if it doesn't exist, insert it
	else if (i->isReal()) {	i->linIx=iOld->linIx; // else replace existing one (with special care to linIx, only valid locally)
		(*this)[i->linIx] = i;}
	return res;
	
// 	assert(bodies);
// 	Body::id_t id1 = i->getId1();
// 	Body::id_t id2 = i->getId2();
// 	if (id1 > id2 ) {swap(id1, id2); } //FIXME: is this really safe? I think it breaks Shop::createExplicitInteraction (Bruno)
// 	
// 	if (!(*bodies)[id1] || !(*bodies)[id2]) {
// 	  LOG_ERROR("Bodies not found in body container --> " << id1 << "  " << id2);
// 	  return false;
// 	}
// 	
// 	const shared_ptr<Body>& b1=(*bodies)[id1];
// 	const shared_ptr<Body>& b2=(*bodies)[id2]; 
// 	
// // 	if (b1 -> getIsSubdomain() || b2 -> getIsSubdomain() ) {return false; }
// 
// 	if(!b1->intrs.insert(Body::MapId2IntrT::value_type(id2,i)).second) return false; // already exists
// 	if(!b2->intrs.insert(Body::MapId2IntrT::value_type(id1,i)).second) return false;
// 	
// 	//check if this interaction already exists : 
// 	linIntrs.resize(++currSize); // currSize updated
// 	linIntrs[currSize-1]=i; // assign last element
// 	i->linIx=currSize-1; // store the index back-reference in the interaction (so that it knows how to erase/move itself
// 	const shared_ptr<Scene>& scene=Omega::instance().getScene();
// 	i->iterMadeReal=scene->iter; 
// 	
// 	
// 
// // 	std::cout << "Interaction contaier size now = " << linIntrs.size() << std::endl; 
// 	return true; 
	
	
	
}



void InteractionContainer::clear(){
	assert(bodies);
	boost::mutex::scoped_lock lock(drawloopmutex);
	FOREACH(const shared_ptr<Body>& b, *bodies) {
		if (b) b->intrs.clear();
	}
	linIntrs.clear();
	currSize=0;
	dirty=true;
}

bool InteractionContainer::erase(Body::id_t id1,Body::id_t id2, int linPos){
	assert(bodies);
	boost::mutex::scoped_lock lock(drawloopmutex);
	if (id1>id2) swap(id1,id2);
	if(id2>=(Body::id_t)bodies->size()) return false;

	const shared_ptr<Body>& b1((*bodies)[id1]);
	const shared_ptr<Body>& b2((*bodies)[id2]);
	int linIx=-1;
	if(!b1 and !b2) {linIx=linPos; LOG_WARN("!b1 and !b2 for id1="+boost::lexical_cast<string>(id1)+" id2="+boost::lexical_cast<string>(id2));}	
	else {
		if (b1) {
			Body::MapId2IntrT::iterator I(b1->intrs.find(id2));
 			if(I==b1->intrs.end()) {LOG_ERROR("I not found for "+boost::lexical_cast<string>(id1)+" "+boost::lexical_cast<string>(id2));
			} else {
				linIx=I->second->linIx;
				assert(linIx==linPos);
				if(linIx!=linPos) { LOG_ERROR("linIx!=linPos"+boost::lexical_cast<string>(id1)+" "+boost::lexical_cast<string>(id2)+" "+boost::lexical_cast<string>(linIx)+" "+boost::lexical_cast<string>(linPos));}
				//erase from body, we also erase from linIntrs below
				b1->intrs.erase(I);
			}
		}
		if (b2) {
			Body::MapId2IntrT::iterator I(b2->intrs.find(id1));
			if(I==b1->intrs.end()) {LOG_ERROR("I not found(2) for "+boost::lexical_cast<string>(id1)+" "+boost::lexical_cast<string>(id2));
			} else {
				linIx=I->second->linIx;
				assert(linIx==linPos);
				if(linIx!=linPos) { LOG_ERROR("linIx!=linPos(2)"+boost::lexical_cast<string>(id1)+" "+boost::lexical_cast<string>(id2)+" "+boost::lexical_cast<string>(linIx)+" "+boost::lexical_cast<string>(linPos));}
				//erase from body, we also erase from linIntrs below
				b2->intrs.erase(I);
			}
		}
	}
	if(linIx<0) {
		LOG_ERROR("InteractionContainer::erase: attempt to delete interaction with a deleted body (the definition of linPos in the call to erase() should fix the problem) for  ##"+boost::lexical_cast<string>(id1)+"+"+boost::lexical_cast<string>(id2));
		return false;}
	// iid is not the last element; we have to move last one to its place
	if (linIx<(int)currSize-1) {
		linIntrs[linIx]=linIntrs[currSize-1];
		linIntrs[linIx]->linIx=linIx; // update the back-reference inside the interaction
	}
	// in either case, last element can be removed now
	linIntrs.resize(--currSize); // currSize updated
	return true;
}

const shared_ptr<Interaction>& InteractionContainer::find(Body::id_t id1,Body::id_t id2){
	assert(bodies);
	if (id1>id2) swap(id1,id2);
	// those checks could be perhaps asserts, but pyInteractionContainer has no access to the body container...
	if(id2>=(Body::id_t)bodies->size()){ empty=shared_ptr<Interaction>(); return empty; }
	const shared_ptr<Body>& b1((*bodies)[id1]);
	if(!b1) { empty=shared_ptr<Interaction>(); return empty; }
	Body::MapId2IntrT::iterator I(b1->intrs.find(id2));
	if (I!=b1->intrs.end()) return I->second;
	else { empty=shared_ptr<Interaction>(); return empty; }
}

bool InteractionContainer::insert(Body::id_t id1,Body::id_t id2)
{
	shared_ptr<Interaction> i(new Interaction(id1,id2) );
	return insert(i);
}

void InteractionContainer::requestErase(Body::id_t id1, Body::id_t id2){
	const shared_ptr<Interaction> I=find(id1,id2); if(!I) return;
	I->reset();
}

void InteractionContainer::requestErase(const shared_ptr<Interaction>& I){
	I->reset();
}

void InteractionContainer::requestErase(Interaction* I){
	I->reset();
}

void InteractionContainer::eraseNonReal(){
	FOREACH(const shared_ptr<Interaction>& i, *this) if(!i->isReal()) this->erase(i->getId1(),i->getId2(),i->linIx);
}

// compare interaction based on their first id
struct compPtrInteraction{
	bool operator() (const shared_ptr<Interaction>& i1, const shared_ptr<Interaction>& i2) const {
		return (*i1)<(*i2);
	}
};

void InteractionContainer::updateSortedIntrs(){
	this->sortedIntrs.resize(this->linIntrs.size());
	this->sortedIntrs = this->linIntrs;
	sort(this->sortedIntrs.begin(), this->sortedIntrs.end(), this->compareTwoInteractions);
};

void InteractionContainer::preSave(InteractionContainer&){
	for(const shared_ptr<Interaction>& I : *this) {
		if(I->geom || I->phys) interaction.push_back(I);
		// since requestErase'd interactions have no interaction physics/geom, they are not saved
	}
	if(serializeSorted) std::sort(interaction.begin(),interaction.end(),compPtrInteraction());
}
void InteractionContainer::postSave(InteractionContainer&){ interaction.clear(); }

void InteractionContainer::preLoad(InteractionContainer&){ interaction.clear(); }

void InteractionContainer::postLoad__calledFromScene(const shared_ptr<BodyContainer>& bb){
	bodies=&bb->body;
	clear();
	FOREACH(const shared_ptr<Interaction>& I, interaction){
		Body::id_t id1=I->getId1(), id2=I->getId2();
		if (!(*bodies)[id1] || !(*bodies)[id2]) {
			return;
		} else {
			insert(I);
		}
	}
	interaction.clear();
}
