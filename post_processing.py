import requests, sys, json
from pickle import load, dump


GO_link = "https://www.ebi.ac.uk/QuickGO/services/ontology/go/terms/"
#criteria = "/ancestors?relations=is_a%2Cpart_of%2Coccurs_in%2Cregulates"
criteria = "/ancestors?relations=is_a"

anchestorList=load(open("GO_ancestorList.p", 'rb'))

# fetch the parents for a GO term
def fetchGoParents(id):
    _,id=id.split(':')
    GOId="GO%3A"+id

    requestURL =GO_link+GOId+criteria

    r = requests.get(requestURL, headers={"Accept": "application/json"})

    if not r.ok:
        return False
    responseBody = r.text
    data=json.loads(responseBody)
    #print data
    if data["results"]==[]:
        return ["GO:"+id]
    isObsolete=data["results"][0]["isObsolete"]
    if isObsolete:
        print (isObsolete)
        return ["GO:"+id]

    parents=data["results"][0]["ancestors"]
    if parents!=[]:
        P=[str(go) for go in parents]

    return P

#iteratively form a list of parent propagating through hierarchy
def dfs_collection(id):
    visited={id:1}
    Ids={id}
    #print id
    while Ids!=set():
        nextid=Ids.pop()

        if visited[nextid]==0:
            if Ids.__contains__(nextid):
                Ids.remove(nextid)
            continue

        Gos=fetchGoParents(nextid)
        if Gos.__contains__(nextid):
            Gos.remove(nextid)
        Ids=set(list(Ids)+Gos)
        #print nextid, Ids

        visited[nextid]=0
        for x in  Gos:
            if x in visited:
                continue
            visited.update({x:1})
        #print visited

    return visited.keys()


def ancestorsViaDAG(predictedNew):
    predictedNewWithAncestor=[]

    for i in predictedNew:
        if i in anchestorList:
            result=anchestorList[i]
        else:
            predictedNewWithAncestor.append(i)
            result=dfs_collection(i)
            anchestorList.update({i:result})
        for val in result:
            predictedNewWithAncestor.append(val)
            

    newset=set(predictedNewWithAncestor)
    predictedNewWithAncestor=list(newset)

    return predictedNewWithAncestor





if __name__ == '__main__':
    x=ancestorsViaDAG(['GO:1903901', 'GO:1903900'])
    print(x)


